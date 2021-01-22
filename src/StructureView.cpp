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
#include <h3dsrc/Text.h>
#include "StructureView.h"
#include "Ensemble.h"
#include "Segment.h"

StructureView::StructureView(QWidget *parent) : SlipGL(parent)
{
	_centreSet = false;
	setBackground(1, 1, 1, 1);
	setZFar(2000.);
	setFocusPolicy(Qt::ClickFocus);
	_text = NULL;
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
	                    0, -300, -20);
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

