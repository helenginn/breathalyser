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

#ifndef __sequins__sequenceview__
#define __sequins__sequenceview__

#include <QMainWindow>

class FastaMaster;
class FastaGroup;
class QBoxLayout;
class QLabel;
class Fasta;
class Main;

class SequenceView : public QMainWindow
{
Q_OBJECT
public:
	SequenceView(QWidget *parent, FastaMaster *m);

	void setMain(Main *main)
	{
		_main = main;
	}
	
	void populate(Fasta *f);
	void populate(FastaGroup *g);
private:
	QLabel *addToLayout(QBoxLayout *box, std::string left, std::string right);
	QBoxLayout *setupWindow();

	FastaMaster *_master;
	Fasta *_fasta;
	FastaGroup *_group;
	Main *_main;
};

#endif
