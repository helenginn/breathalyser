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

#ifndef __breathalyser__loadfastas__
#define __breathalyser__loadfastas__

#include <QMainWindow>

class Main;
class QLineEdit;
class QCheckBox;

class LoadFastas : public QMainWindow
{
Q_OBJECT
public:
	LoadFastas(QWidget *parent = NULL);
	void setMain(Main *m);

	void setProtein(bool p);
	void loadFastas(std::string filename, int start, int end);
	void loadSequence(std::string filename, int start, int end, 
	                  bool isProtein);
public slots:
	void loadChosenFasta();
	void chooseFasta();
private:
	QLineEdit *_fastaLine;
	QLineEdit *_rangeLine;
	QCheckBox *_isProtein;
	
	Main *_main;
};

#endif
