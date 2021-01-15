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

#include <QMenu>
#include <QMenuBar>
#include <iostream>

#include "Main.h"
#include "LoadStructure.h"
#include "Ensemble.h"
#include "StructureView.h"

Main::Main(QWidget *parent) : QMainWindow(parent)
{
	setGeometry(0, 0, 1500, 900);
	_pdbTree = new QTreeWidget(this);
	_pdbTree->show();
	_pdbTree->setHeaderLabel("Structures");
	
	connect(_pdbTree, &QTreeWidget::currentItemChanged,
	        this, &Main::clickedStructure);

	_diffTree = new QTreeWidget(this);
	_diffTree->show();
	_diffTree->setHeaderLabel("Distance matrices");
	_view = new StructureView(this);
	_view->show();
	makeMenu();
	show();
}

void Main::makeMenu()
{
// QList<QPushButton *> allPButtons = parentWidget.findChildren<QPushButton *>();

	QMenu *structures = menuBar()->addMenu(tr("&Structures"));
	QAction *act = structures->addAction(tr("&Load structures"));
	connect(act, &QAction::triggered, this, &Main::loadStructures);
}

void Main::resizeEvent(QResizeEvent *e)
{
	int x = 0;
	int y = menuBar()->height();
	int w = width();
	int h = height() - y;

	_pdbTree->setGeometry(x, y, 250, h / 2);
	_diffTree->setGeometry(x, y + h / 2, 250, h / 2);
	_view->setGeometry(x + 250, y, w - 500, h - 100);
}

void Main::loadStructures()
{
	LoadStructure *load = new LoadStructure(NULL);
	load->setMain(this);
}

void Main::receiveEnsemble(Ensemble *e)
{
	_pdbTree->addTopLevelItem(e);
	
}

void Main::clickedStructure()
{
	QTreeWidgetItem *item = _pdbTree->currentItem();
	Ensemble *e = dynamic_cast<Ensemble *>(item);
	
	if (!e)
	{
		std::cout << "Not an ensemble!" << std::endl;
		return;
	}
	
	_view->addEnsemble(e);
}
