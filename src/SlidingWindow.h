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

#ifndef __breathalyser__slidingwindow__
#define __breathalyser__slidingwindow__

#include <QMainWindow>

class QCheckBox;
class QLineEdit;
class QSlider;
class QLabel;
class Main;

class SlidingWindow : public QMainWindow
{
Q_OBJECT
public:
	SlidingWindow(QWidget *parent, Main *main);

	void setMain(Main *m)
	{
		_main = m;
	}

public slots:
	void run();
	void slid(int val);
private:
	QLineEdit *_mutLine;
	QSlider *_genomeCount;
	QLabel *_genomeReport;
	QCheckBox *_overwrite;
	int _scale;

	Main *_main;
};

#endif
