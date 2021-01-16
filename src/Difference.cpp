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

#include "Ensemble.h"
#include "Difference.h"
#include <QPainter>
#include <libsrc/Polymer.h>
#include <libsrc/Atom.h>

Difference::Difference(int w, int h) : QImage(w + 1, h, QImage::Format_RGB32)
{

}

void Difference::findCommonAtoms()
{
	if (_ea == NULL || _eb == NULL || 
	    _ea->crystal() == NULL ||
	    _eb->crystal() == NULL)
	{
		return;
	}

	for (size_t i = 0; i < _ea->chainCount(); i++)
	{
		std::string ch = _ea->chain(i);
		
		AtomList as = _ea->crystal()->findAtoms("CA", INT_MAX, ch);
		
		for (size_t j = 0; j < as.size(); j++)
		{
			int resNum = as[j]->getResidueNum();
			AtomList bs = _eb->crystal()->findAtoms("CA", resNum, ch);
			AtomCouple couple;
			
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

void Difference::populate()
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
		double l1 = NAN;
		if (c1.second)
		{
			vec3 a = c1.first->getInitialPosition();
			vec3 b = c1.second->getInitialPosition();
			vec3 diff = vec3_subtract_vec3(b, a);
			l1 = vec3_length(diff);
		}

		for (int i = 0; i < num; i++)
		{
			AtomCouple c2 = _atomCouples[i];
			double l2 = NAN;
			if (c2.second)
			{
				vec3 a = c2.first->getInitialPosition();
				vec3 b = c2.second->getInitialPosition();
				vec3 diff = vec3_subtract_vec3(b, a);
				l2 = vec3_length(diff);
			}
			
			double val = l2 - l1;
			
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
}
