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

#ifndef __breathalyser__diffdisplay__
#define __breathalyser__diffdisplay__

#include <QWidget>

class Main;
class QLabel;
class QSlider;
class QPushButton;
class Difference;

class DiffDisplay : public QWidget
{
public:
	DiffDisplay(QWidget *parent, Difference *diff);

	void changeDifference(Difference *diff = NULL);
	
	Difference *difference()
	{
		return _diff;
	}
	
	void setMain(Main *main)
	{
		_main = main;
	}
private:
	Main *_main;
	Difference *_diff;
	QLabel *_label;
	QSlider *_threshold;
	QPushButton *_calculate;
};

#endif
