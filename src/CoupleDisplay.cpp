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

#include <QVBoxLayout>
#include <QPushButton>
#include "CoupleDisplay.h"
#include "StructureView.h"
#include "Difference.h"
#include "Ensemble.h"

CoupleDisplay::CoupleDisplay(QWidget *parent, StructureView *view)
: QWidget(parent)
{
	_view = view;
	_type = DisplayBoth;
	_main = NULL;
	view->setParent(this);

	QVBoxLayout *box = new QVBoxLayout();
	box->addWidget(_view);

	QHBoxLayout *hbox = new QHBoxLayout();

	QPushButton *f = new QPushButton("First structure", this);
	QPushButton *s = new QPushButton("Second structure", this);
	QPushButton *b = new QPushButton("Both structures", this);
	QPushButton *m = new QPushButton("Merge segments", this);
	hbox->addWidget(f);
	hbox->addWidget(s);
	hbox->addWidget(b);
	hbox->addWidget(m);
	
	connect(f, &QPushButton::clicked, this, &CoupleDisplay::displayFirst);
	connect(s, &QPushButton::clicked, this, &CoupleDisplay::displaySecond);
	connect(b, &QPushButton::clicked, this, &CoupleDisplay::displayBoth);
	connect(m, &QPushButton::clicked, this, &CoupleDisplay::mergeSegments);

	box->addLayout(hbox);
	setLayout(box);
	
	_ea = NULL;
	_eb = NULL;

}

void CoupleDisplay::setDifference(Difference *diff)
{
	diff->applySegmentsToEnsembles();
	_diff = diff;
	setEnsembles(diff->aEnsemble(), diff->bEnsemble());
}

void CoupleDisplay::setEnsembles(Ensemble *a, Ensemble *b)
{
	_view->clearObjects();

	_ea = a;
	_eb = b;

	_view->addEnsemble(_ea);
	_view->addEnsemble(_eb);

	correctDisplay();
}

void CoupleDisplay::correctDisplay()
{
	switch (_type)
	{
		case DisplayFirst:
		displayFirst();
		break;

		case DisplaySecond:
		displaySecond();
		break;

		case DisplayBoth:
		displayBoth();
		break;

		default:
		break;
	}
}

void CoupleDisplay::displayFirst()
{
	_view->clearObjects();
	_view->addEnsemble(_ea);
	_type = DisplayFirst;
}

void CoupleDisplay::displaySecond()
{
	_view->clearObjects();
	_view->addEnsemble(_eb);
	_type = DisplaySecond;
}

void CoupleDisplay::displayBoth()
{
	_view->clearObjects();
	_view->addEnsemble(_ea);
	_view->addEnsemble(_eb);
	_type = DisplayBoth;
}

void CoupleDisplay::mergeSegments()
{
	_diff->tryMerges();

}
