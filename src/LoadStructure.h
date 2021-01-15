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

#ifndef __breathalyser__loadstructure__
#define __breathalyser__loadstructure__

#include <QMainWindow>

class QLineEdit;
class Main;

class LoadStructure : public QMainWindow
{
Q_OBJECT
public:
	LoadStructure(QWidget *parent = NULL);

	void setMain(Main *m)
	{
		_main = m;
	}

public slots:
	void choosePDB();
	void loadPDB();
private:
	QLineEdit *_pdbLine;
	Main *_main;

};

#endif
