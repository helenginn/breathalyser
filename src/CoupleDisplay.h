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

#ifndef __breathalyser__coupledisplay__
#define __breathalyser__coupledisplay__

#include <QWidget>

class Main;
class Ensemble;
class Difference;
class StructureView;

class CoupleDisplay : public QWidget
{
Q_OBJECT
public:
	CoupleDisplay(QWidget *parent, StructureView *view);

	void setMain(Main *main)
	{
		_main = main;
	}

	void setDifference(Difference *d);
public slots:
	void displayFirst();
	void displaySecond();
	void displayBoth();
	void mergeSegments();

private:
	void setEnsembles(Ensemble *a, Ensemble *b);
	void correctDisplay();

	typedef enum
	{
		DisplayFirst,
		DisplaySecond,
		DisplayBoth
	} DisplayType;
	
	DisplayType _type;

	Main *_main;
	StructureView *_view;

	Ensemble *_ea;
	Ensemble *_eb;
	Difference *_diff;
};

#endif
