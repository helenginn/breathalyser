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
#include <QTabWidget>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <iostream>

#include "DiffDisplay.h"
#include "Difference.h"
#include "Main.h"
#include "MyDictator.h"
#include "LoadStructure.h"
#include "Ensemble.h"
#include "StructureView.h"

Main::Main(QWidget *parent) : QMainWindow(parent)
{
	makeMenu();
	setGeometry(0, 0, 1300, 900);
	
	_ref = NULL;

	QWidget *window = new QWidget();
	QHBoxLayout *layout = new QHBoxLayout();
	window->setLayout(layout);

	QVBoxLayout *treeout = new QVBoxLayout();

	_pdbTree = new QTreeWidget(window);
	_pdbTree->show();
	_pdbTree->setHeaderLabel("Structures");
	_pdbTree->setSelectionMode(QAbstractItemView::ExtendedSelection);
	_pdbTree->setContextMenuPolicy(Qt::CustomContextMenu);
	_pdbTree->setMaximumSize(QSize(250, 500));

	_diffTree = new QTreeWidget(window);
	_diffTree->show();
	_diffTree->setHeaderLabel("Distance matrices");
	_diffTree->setMaximumSize(QSize(250, 500));

	treeout->addWidget(_pdbTree);
	treeout->addWidget(_diffTree);
	layout->addItem(treeout);

	_tabs = new QTabWidget(window);
	layout->addWidget(_tabs);
	
	setCentralWidget(window);

	_view = new StructureView(NULL);
	_tabs->addTab(_view, "Structure view");

	_diff = new DiffDisplay(NULL, NULL);
	_diff->setMain(this);
	_tabs->addTab(_diff, "Difference view");

	_couple = new StructureView(NULL);
	_tabs->addTab(_couple, "Couple view");

	connect(_pdbTree, &QTreeWidget::customContextMenuRequested,
	        this, &Main::structureMenu);
	
	connect(_pdbTree, &QTreeWidget::itemSelectionChanged,
	        this, &Main::clickedStructure);
	
	connect(_diffTree, &QTreeWidget::currentItemChanged,
	        this, &Main::clickedDifference);

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
	_diff->changeDifference();
}

void Main::loadStructures()
{
	LoadStructure *load = new LoadStructure(NULL);
	load->setMain(this);
	load->show();
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
	
	d->toCoupleView(_couple);
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
	diff->setMain(this);
	diff->setEnsembles(e1, e2);
	diff->populate();
	_diff->changeDifference(diff);
	_diffTree->addTopLevelItem(diff);
	_diffTree->setCurrentItem(diff);
	
	_tabs->setCurrentWidget(_diff);
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

void Main::setCommandLineArgs(int argc, char *argv[])
{
	for (int i = 1; i < argc; i++)
	{
		std::string str = argv[i];
		_args.push_back(str);
	}

	_dictator = new MyDictator(this);
	_dictator->setArgs(_args);
	_dictator->run();
}

