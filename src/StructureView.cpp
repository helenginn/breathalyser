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

#include <iostream>
#include <QMenu>
#include <h3dsrc/Text.h>
#include <hcsrc/FileReader.h>
#include "StructureView.h"
#include "Ensemble.h"
#include "Segment.h"
#include "FastaMaster.h"
#include "Main.h"

StructureView::StructureView(QWidget *parent) : SlipGL(parent)
{
	_centreSet = false;
	setBackground(1, 1, 1, 1);
	setZFar(2000.);
	setFocusPolicy(Qt::ClickFocus);
	_text = NULL;
	_ensemble = NULL;
}

void StructureView::initializeGL()
{
	SlipGL::initializeGL();
}

void StructureView::addEnsemble(Ensemble *e)
{
	if (e == NULL)
	{
		return;
	}

	e->repopulate();

	addObject(e, false);
	_ensemble = e;
	
	if (!_centreSet)
	{
		vec3 centre = e->averagePos();
		_centre = centre;
		focusOnPosition(centre);
		_centreSet = true;
	}
}

void StructureView::addLabel(std::string str)
{
	Text *text = new Text();
	text->setProperties(_centre, str, 72, Qt::black,
	                    0, -33, -20);
	text->prepare();
	setText(text);
}

void StructureView::setText(Text *text)
{
	removeObject(_text);
	delete _text;
	_text = text;
	if (text)
	{
		addObject(text, false);
	}
}

void StructureView::clickMouse(double x, double y)
{
	if (_ensemble != NULL)
	{
		std::string mut = _ensemble->whichMutation(x, y);
		std::cout << "Selected residue: " << mut << std::endl;
	}
}

void StructureView::mouseReleaseEvent(QMouseEvent *e)
{
	if (!_moving && e->button() == Qt::LeftButton)
	{
		double x = e->x(); double y = e->y();
		convertCoords(&x, &y);
		clickMouse(x, y);
	}
	if (!_moving && e->button() == Qt::RightButton)
	{
		double x = e->x(); double y = e->y();
		QPoint p = mapToGlobal(QPoint(x, y));
		makeMutationMenu(p);
	}

	SlipGL::mouseReleaseEvent(e);
}

void StructureView::makeMutationMenu(QPoint &p)
{
	std::string resi = _ensemble->selectedMutation();
	std::string inverse = "!" + resi;

	if (resi.length() == 0)
	{
		return;
	}
	
	QMenu *m = new QMenu();
	QAction *act = m->addAction("Select on mutation");
	connect(act, &QAction::triggered, 
	        this, [=] {_main->fMaster()->requireMutation(resi);});
	act = m->addAction("Select on wild-type");
	connect(act, &QAction::triggered, 
	        this, [=] {_main->fMaster()->requireMutation(inverse);});
	m->exec(p);
}

void StructureView::screenshot(std::string filename)
{
	std::string zero = getBaseFilenameWithPath(filename) + "_0.";
	zero += getExtension(filename);

	std::string one = getBaseFilenameWithPath(filename) + "_1.";
	one += getExtension(filename);

	std::cout << zero << std::endl;
	std::cout << one << std::endl;

	_ensemble->setMode(0);
	update();
	saveImage(zero);
	_ensemble->setMode(1);
	update();
	saveImage(one);
	_ensemble->setMode(-1);
	update();
	saveImage(filename);
	std::cout << filename << std::endl;

}
