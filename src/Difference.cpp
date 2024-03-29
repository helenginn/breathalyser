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
#include "CoupleDisplay.h"
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
	_coupleDisplay = NULL;
}

void Difference::toCoupleDisplay(CoupleDisplay *view)
{
	_coupleDisplay = view;
	_coupleDisplay->setDifference(this);
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
				_map[as[j]] = bs[0];
				_map[bs[0]] = as[j];
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
			
			vec3 diff = vec3_subtract_vec3(b1, b2);
			double l1 = vec3_length(diff);
			diff = vec3_subtract_vec3(a1, a2);
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

	QColor c = QColor(0, 0, 0, 255);
	QPen p = QPen(c);
	p.setWidth(box_size);
	QBrush b = QBrush(c, Qt::BDiagPattern);
	painter.setPen(p);
	painter.setBrush(b);
	
	applySegmentsToEnsembles();
	
	for (size_t i = 0; i < _ea->segmentCount(); i++)
	{
		Segment *s = _ea->segment(i);
		
		for (size_t j = 0; j < s->subCount(); j++)
		{
			for (size_t k = 0; k < s->subCount(); k++)
			{
				drawOnSegment(s->sub(j), s->sub(k), painter, box_size);
			}
		}
	}
	
	_drawn = true;
}

void Difference::drawOnSegment(Segment *s, Segment *t, 
                               QPainter &painter, double box_size)
{

	AtomPtr first = _atomCouples[0].first;
	std::string ch = first->getChainID().substr(0, 1);
	int min, max;
	_ea->minMaxResidues(ch, &min, &max);

	int ss, se;
	s->startEnd(&ss, &se);
	ss -= min;
	se -= min;
	int ssize = se - ss;
	
	painter.drawRect(box_size * ss, box_size * ss,
	                 box_size * ssize, box_size * ssize);
	
	if (s == t)
	{
		return;
	}

	int ts, te;
	t->startEnd(&ts, &te);
	ts -= min;
	te -= min;
	int tsize = te - ts;
	
	painter.drawRect(box_size * ss, box_size * ts,
	                 box_size * ssize, box_size * tsize);

}

void Difference::addSegments(Segment *seg_a, double threshold)
{
	Segment *seg_b = new Segment(this, threshold);

	for (size_t i = 0; i < seg_a->atomCount(); i++)
	{
		AtomPtr a = seg_a->atom(i);
		AtomPtr b = _map[a];
		
		if (b)
		{
			seg_b->addAtom(b);
		}
	}

	seg_a->setSister(seg_b);
	seg_b->setSister(seg_a);

	_ea->addSegment(seg_a);
	_eb->addSegment(seg_b);
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
	v *= 5.0;

	double real = v / _max;

	_ea->deleteSegments();
	_eb->deleteSegments();

	Segment *seg_a = new Segment(this, real);

	AtomPtr start_atom;
	for (size_t i = 0; i < _atoms.size(); i++)
	{
		AtomPtr a = _atoms[i];
		
		bool ok = seg_a->addAtomIfValid(a);
		
		if (!ok || i == _atoms.size() - 1)
		{
			if (seg_a->atomCount() >= 3)
			{
				addSegments(seg_a, real);
				seg_a = new Segment(this, real);
			}
			else
			{
				delete seg_a;
				seg_a = new Segment(this, real);
			}
		}
	}

	_aSegments = _ea->segments();
	_bSegments = _eb->segments();
}

void Difference::applySegmentsToEnsembles()
{
	_ea->setSegments(_aSegments);
	_eb->setSegments(_bSegments);

}

void Difference::thresholdChanged(int val)
{
	findSegments(val);
	populate();
	_display->changeDifference(NULL);
}

void Difference::calculate()
{
	for (size_t i = 0; i < _ea->segmentCount(); i++)
	{
		_ea->segment(i)->populate();
		_eb->segment(i)->populate();
		_ea->segment(i)->forceSisterColour();
		_ea->segment(i)->refineMesh();
		_eb->segment(i)->refineMesh();
	}
}

double Difference::valueBetween(AtomPtr a, AtomPtr b)
{
	AtomCouple c = std::make_pair(a, b);
	if (!_vals.count(c))
	{
		return NAN;
	}
	
	return _vals[c];
}

bool Difference::closerSegmentThan(Segment *s, Segment *t)
{
	vec3 as = s->AtomGroup::centroid();
	vec3 at = t->AtomGroup::centroid();
	vec3 targ_diff = vec3_subtract_vec3(as, at);
	double tl = vec3_length(targ_diff);

	for (size_t i = 0; i < _ea->segmentCount(); i++)
	{
		Segment *a = _ea->segment(i);

		if (a == s || a == t)
		{
			continue;
		}
		
		vec3 ac = a->AtomGroup::centroid();
		vec3 t1 = vec3_subtract_vec3(ac, as);
		vec3 t2 = vec3_subtract_vec3(ac, at);
		
		if (vec3_length(t1) < tl || vec3_length(t2) < tl)
		{
			return true;
		}
	}
	
	return false;
}

void Difference::tryMerges()
{
	double mult = 2.0;
	applySegmentsToEnsembles();
	for (size_t i = 0; i < _ea->segmentCount() - 1; i++)
	{
		Segment *a1 = _ea->segment(i);

		for (size_t j = i + 1; j < _ea->segmentCount(); j++)
		{
			Segment *a2 = _ea->segment(j);

			if (closerSegmentThan(a1, a2))
			{
				continue;
			}

			if (!a1->enoughCommonGround(a2, mult))
			{
				continue;
			}

			std::cout << "Combined two segments " << i << 
			" and " << j << std::endl;

			Segment *b1 = _eb->segment(i);
			Segment *b2 = _eb->segment(j);

			Segment *aMaster = Segment::segmentFrom(a1, a2, mult);
			Segment *bMaster = Segment::segmentFrom(b1, b2, mult);

			aMaster->setSister(bMaster);
			bMaster->setSister(aMaster);

			_ea->removeSegment(i);
			_eb->removeSegment(i);
			i--;
			j--;
			_ea->removeSegment(j);
			_eb->removeSegment(j);
			j--;

			_ea->addSegment(aMaster);
			_eb->addSegment(bMaster);
			break;
		}
	}

	_aSegments = _ea->segments();
	_bSegments = _eb->segments();

	calculate();
	populate();
	_display->changeDifference(NULL);
}

AtomPtr Difference::getPartnerAtom(AtomPtr a)
{
	return _map[a];
}
