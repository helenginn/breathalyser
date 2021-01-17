// breathalyser
// Copyright (C) 2019 Helen Ginn
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
// 
// Please email: vagabond @ hginn.co.uk for more details.

#include "Segment.h"
#include "Ensemble.h"
#include "Difference.h"
#include "DiffDisplay.h"
#include "StructureView.h"
#include "Main.h"
#include <QPainter>
#include <libsrc/Polymer.h>
#include <libsrc/Atom.h>
#include <iomanip>

Difference::Difference(int w, int h) : QImage(w + 1, h, QImage::Format_RGB32)
{
	_drawn = false;
	_display = NULL;
	_main = NULL;
	_max = 0;
	_ea = NULL;
	_eb = NULL;
	_coupleView = NULL;
}

void Difference::toCoupleView(StructureView *view)
{
	_coupleView = view;
	view->clearObjects();
	
	view->addEnsemble(_ea);
	view->addEnsemble(_eb);
}

void Difference::findCommonAtoms()
{
	if (_ea == NULL || _eb == NULL || 
	    _ea->crystal() == NULL ||
	    _eb->crystal() == NULL ||
	    _main == NULL)
	{
		return;
	}
	
	Ensemble *ref = _main->reference();
	_atoms.clear();

	for (size_t i = 0; i < ref->chainCount(); i++)
	{
		std::string ch = ref->chain(i);
		std::string cha = ref->findMatchingChain(ch, _ea);
		std::string chb = ref->findMatchingChain(ch, _eb);
		
		AtomList as = _ea->crystal()->findAtoms("CA", INT_MAX, cha);
		
		for (size_t j = 0; j < as.size(); j++)
		{
			int resNum = as[j]->getResidueNum();
			AtomList bs = _eb->crystal()->findAtoms("CA", resNum, chb);
			AtomCouple couple;
			_atoms.push_back(as[j]);
			
			if (bs.size() >= 1)
			{
				couple = std::make_pair(as[j], bs[0]);
			}
			else
			{
				couple = std::make_pair(as[j], AtomPtr());
			}
			
			_atomCouples.push_back(couple);
		}
	}
}

void Difference::setEnsembles(Ensemble *a, Ensemble *b)
{
	_ea = a;
	_eb = b;
	
	std::string name = a->name() + " to " + b->name();
	this->QTreeWidgetItem::setText(0, QString::fromStdString(name));
	
	findCommonAtoms();
}

void Difference::populate(bool force)
{
	int num = _atomCouples.size();

	QPainter painter(this);

	double box_size = ((double)width() / (double)(num));
	
	int red = 255;
	int green = 0;
	int blue = 0;

	std::string chain = "";
	for (int j = 0; j < num; j++)
	{
		AtomCouple c1 = _atomCouples[j];
		vec3 a1 = make_vec3(NAN, NAN, NAN);
		vec3 b1 = make_vec3(NAN, NAN, NAN);

		if (c1.second)
		{
			a1 = c1.first->getInitialPosition();
			b1 = c1.second->getInitialPosition();
		}

		for (int i = 0; i < num; i++)
		{
			AtomCouple c2 = _atomCouples[i];
			vec3 a2 = make_vec3(NAN, NAN, NAN);
			vec3 b2 = make_vec3(NAN, NAN, NAN);

			if (c2.second)
			{
				a2 = c2.first->getInitialPosition();
				b2 = c2.second->getInitialPosition();
			}
			
			vec3 diff = vec3_subtract_vec3(b1, a1);
			double l1 = vec3_length(diff);
			diff = vec3_subtract_vec3(b2, a2);
			double l2 = vec3_length(diff);

			double val = l2 - l1;

			if (!_drawn)
			{
				AtomCouple x1 = std::make_pair(c1.first, c2.first);
				AtomCouple x2 = std::make_pair(c1.second, c2.second);

				_vals[x1] = val;
				_vals[x2] = val;
			}
			
			if (fabs(val) > _max)
			{
				_max = fabs(val);
			}
			
			if (val > 2) val = 2;
			if (val < -2) val = -2;

			if (val != val) /* we go grey */
			{
				red = 100;
				green = 100;
				blue = 100;
			}
			else if (val <= -1) /* we go black */
			{
				val = -(val + 1.);
				red = 0;
				green = 0;
				blue = 255 - val * 255;
			}
			else if (val < 0)
			{
				/* we go blue. */
				val = -val;
				red = 255 - val * 255;
				green = 255 - val * 255;
				blue = 255;
			}
			else if (val >= 1.0) /* We go yellow. */
			{
				val -= 1; 
				red = 255;
				green = val * 255;
				blue = 0;
			}
			else if (val >= 0) /* We go red. */
			{
				red = 255;
				green = 255 - val * 255;
				blue = 255 - val * 255;
			}

			QColor c = QColor(red, green, blue, 255);
			QPen p = QPen(c);
			QBrush b = QBrush(c, Qt::SolidPattern);
			painter.setPen(p);
			painter.setBrush(b);
			
			painter.drawRect(box_size * i, box_size * j,
			                 box_size + 1, box_size + 1);
		}
	}
	
	for (size_t i = 0; i < _segments.size(); i++)
	{
		QColor c = QColor(0, 0, 0, 255);
		QPen p = QPen(c);
		p.setWidth(box_size);
		QBrush b = QBrush(c, Qt::BDiagPattern);
		painter.setPen(p);
		painter.setBrush(b);
		
		AtomPtr first = _atomCouples[0].first;
		std::string ch = first->getChainID().substr(0, 1);
		int min, max;
		_ea->minMaxResidues(ch, &min, &max);

		int start, end;
		_segments[i]->startEnd(&start, &end);
		start -= min;
		end -= min;
		int size = end - start;

		painter.drawRect(box_size * start, box_size * start,
		                 box_size * size, box_size * size);
	}
	
	_drawn = true;
}

void Difference::findSegments(double val)
{
	if (val == 0)
	{
		return;
	}

	QSlider *s = static_cast<QSlider *>(QObject::sender());

	double v = val;
	double max = s->maximum();
	v = v / max;
	v = (1 - v);
	v *= 5.0;

	double real = v / _max;

	Segment *seg = new Segment();
	for (size_t i = 0; i < _segments.size(); i++)
	{
		delete _segments[i];
	}

	_segments.clear();

	AtomPtr start_atom;
	for (size_t i = 0; i < _atoms.size(); i++)
	{
		AtomPtr a = _atoms[i];
		bool ok = true;
		
		if (!a)
		{
			ok = false;
		}
		else if (a && seg->atomCount() == 0)
		{
			seg->addAtom(a);
			start_atom = a;
		}
		else if (a && seg->atomCount() >= 1)
		{
			AtomPtr last = _atoms[i-1];
			
			if (last->getResidueNum() != a->getResidueNum() - 1)
			{
				ok = false;
			}

			for (size_t j = 0; j < seg->atomCount(); j++)
			{
				AtomPtr other = seg->atom(j);
				AtomCouple c = std::make_pair(other, a);
				if (!_vals.count(c))
				{
					ok = false;
					break;
				}
				else
				{
					double val = _vals[c];
					if (val > real || val != val)
					{
						ok = false;
						break;
					}
				}
			}
		}
		
		if (ok == false || i == _atoms.size() - 1)
		{
			if (seg->atomCount() > 3)
			{
				_segments.push_back(seg);
				seg = new Segment();
			}
			else
			{
				delete seg;
				seg = new Segment();
			}
		}
		else if (ok)
		{
			seg->addAtom(a);
		}
	}

}

void Difference::thresholdChanged(int val)
{
	findSegments(val);
	populate();
	_display->changeDifference(NULL);
}

void Difference::loadSegmentsToView()
{
	if (!_coupleView)
	{
		return;
	}

	_coupleView->clearSegments();

	for (size_t i = 0; i < segmentCount(); i++)
	{
		_coupleView->addObject(segment(i), false);
	}
}

void Difference::calculate()
{
	loadSegmentsToView();

	for (size_t i = 0; i < segmentCount(); i++)
	{
		segment(i)->populate();
	}
}
