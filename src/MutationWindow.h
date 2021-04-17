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

#ifndef __breathalyser__mutationwindow__
#define __breathalyser__mutationwindow__

#include <QMainWindow>

class QLineEdit;
class Main;

class MutationWindow : public QMainWindow
{
Q_OBJECT
public:
	MutationWindow(QWidget *parent);

	void mutationWindow();

	void setMain(Main *m)
	{
		_main = m;
	}

	void setScan(bool sc);

public slots:
	void run();
private:
	Main *_main;
	QLineEdit *_mutLine;
	bool _scan;

};

#endif
