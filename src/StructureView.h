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

#ifndef __breathalyser__structureview__
#define __breathalyser__structureview__

#include <QObject>
#include <h3dsrc/SlipGL.h>

class Main;
class Ensemble;
class Text;

class StructureView : public SlipGL
{
public:
	StructureView(QWidget *parent);
	void addEnsemble(Ensemble *e);
	
	void setText(Text *text);
	void addLabel(std::string string);
	
	void setMain(Main *main)
	{
		_main = main;
	}

	void screenshot(std::string filename);
protected:
	void makeMutationMenu(QPoint &p);
	void clickMouse(double x, double y);

	virtual void initializeGL();
	virtual void mouseReleaseEvent(QMouseEvent *e);
private:
	Main *_main;
	Ensemble *_ensemble;
	Text *_text;
	vec3 _centre;
	bool _centreSet;

};

#endif
