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

#include "DiffDisplay.h"
#include "Difference.h"
#include "Main.h"
#include "LoadStructure.h"
#include "Ensemble.h"
#include "StructureView.h"

Main::Main(QWidget *parent) : QMainWindow(parent)
{
	setGeometry(0, 0, 1300, 900);
	
	_ref = NULL;

	_pdbTree = new QTreeWidget(this);
	_pdbTree->show();
	_pdbTree->setHeaderLabel("Structures");
	_pdbTree->setSelectionMode(QAbstractItemView::ExtendedSelection);
	_pdbTree->setContextMenuPolicy(Qt::CustomContextMenu);
	connect(_pdbTree, &QTreeWidget::customContextMenuRequested,
	        this, &Main::structureMenu);
	
	connect(_pdbTree, &QTreeWidget::itemSelectionChanged,
	        this, &Main::clickedStructure);

	_diffTree = new QTreeWidget(this);
	_diffTree->show();
	_diffTree->setHeaderLabel("Distance matrices");
	
	connect(_diffTree, &QTreeWidget::currentItemChanged,
	        this, &Main::clickedDifference);

	_view = new StructureView(this);
	_view->show();
	_diff = new DiffDisplay(this, NULL);
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
	_diff->setGeometry(x + 250, y, w - 500, h - 100);
	_diff->changeDifference();
}

void Main::loadStructures()
{
	LoadStructure *load = new LoadStructure(NULL);
	load->setMain(this);
}

void Main::receiveEnsemble(Ensemble *e)
{
	_pdbTree->addTopLevelItem(e);
	
	if (_pdbTree->topLevelItemCount() == 1)
	{
		makeReference(e);
	}
}

void Main::clickedStructure()
{
	QList<QTreeWidgetItem *> list = _pdbTree->selectedItems();
	
	_diff->hide();
	_view->show();
	_view->clearObjects();

	for (int i = 0; i < list.size(); i++)
	{
		QTreeWidgetItem *item = list[i];
		Ensemble *e = dynamic_cast<Ensemble *>(item);

		if (!e)
		{
			continue;
		}

		_view->addEnsemble(e);
	}
}

void Main::clickedDifference()
{
	QTreeWidgetItem *item = _diffTree->currentItem();
	Difference *d = dynamic_cast<Difference *>(item);
	
	if (!d)
	{
		std::cout << "Not a difference matrix!" << std::endl;
		return;
	}
	
	_view->hide();
	_diff->show();
	_diff->changeDifference(d);

}

void Main::structureMenu(const QPoint &p)
{
	QMenu *m = new QMenu();
	bool keep = false;

	if (_pdbTree->selectedItems().size() == 2)
	{
		keep = true;
		QAction *act = m->addAction("Difference matrix");
		connect(act, &QAction::triggered, this, &Main::makeDifference);
	}
	
	if (_pdbTree->selectedItems().size() == 1)
	{
		keep = true;
		QAction *act = m->addAction("Set as reference");
		connect(act, &QAction::triggered, this, &Main::setChosenAsReference);
	}
	
	if (keep)
	{
		QPoint tmp = p;
		tmp.setY(tmp.y() + menuBar()->height());
		m->exec(mapToGlobal(tmp));
	}
	else
	{
		delete m;
	}
}

void Main::makeDifference()
{
	QList<QTreeWidgetItem *> list = _pdbTree->selectedItems();
	if (list.size() < 2)
	{
		return;
	}

	Ensemble *e1 = dynamic_cast<Ensemble *>(list[0]);
	Ensemble *e2 = dynamic_cast<Ensemble *>(list[1]);

	if (e1 == NULL || e2 == NULL)
	{
		return;
	}
	
	Difference *diff = new Difference(1000, 1000);
	diff->setEnsembles(e1, e2);
	diff->populate();
	_diffTree->addTopLevelItem(diff);
	_diffTree->setCurrentItem(diff);
}

Ensemble *Main::ensemble(int i)
{
	QTreeWidgetItem *item = _pdbTree->topLevelItem(i);
	return dynamic_cast<Ensemble *>(item);
}

void Main::makeReference(Ensemble *e)
{
	for (int i = 0; i < ensembleCount(); i++)
	{
		ensemble(i)->setReference(false);
	}

	e->setReference(true);
	_ref = e;
}

void Main::setChosenAsReference()
{
	QList<QTreeWidgetItem *> list = _pdbTree->selectedItems();
	if (list.size() < 1)
	{
		return;
	}

	Ensemble *e = dynamic_cast<Ensemble *>(list[0]);
	makeReference(e);
}
