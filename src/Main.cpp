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
#include <fstream>

#include <h3dsrc/Dialogue.h>
#include <h3dsrc/CurveView.h>

#include "SlidingWindow.h"
#include "FastaMaster.h"
#include "FastaGroup.h"
#include "DiffDisplay.h"
#include "CoupleDisplay.h"
#include "Difference.h"
#include "Main.h"
#include "MyDictator.h"
#include "LoadStructure.h"
#include "MutationWindow.h"
#include "SequenceView.h"
#include "LoadFastas.h"
#include "Fasta.h"
#include "Ensemble.h"
#include "StructureView.h"

Main::Main(QWidget *parent) : QMainWindow(parent)
{
	setGeometry(0, 0, 1300, 900);
	
	_ref = NULL;
	_sw = NULL;

	QWidget *window = new QWidget();
	QHBoxLayout *layout = new QHBoxLayout();
	window->setLayout(layout);

	QVBoxLayout *treeout = new QVBoxLayout();

	_pdbTree = new QTreeWidget(NULL);
	_pdbTree->setHeaderLabel("Structures");
	_pdbTree->setSelectionMode(QAbstractItemView::ExtendedSelection);
	_pdbTree->setContextMenuPolicy(Qt::CustomContextMenu);
	_pdbTree->setMaximumSize(QSize(250, 500));

	_diffTree = new QTreeWidget(NULL);
	_diffTree->setHeaderLabel("Distance matrices");
	_diffTree->setMaximumSize(QSize(250, 500));
	
	_fMaster = new FastaMaster(window);
	_fMaster->setHeaderLabel("Sequence groups");
	_fMaster->setMaximumSize(QSize(250, 1500));
	_fMaster->setContextMenuPolicy(Qt::CustomContextMenu);
	_fMaster->setSelectionMode(QAbstractItemView::ExtendedSelection);
	treeout->addWidget(_fMaster);

	/*
	_curveTree = new QTreeWidget(NULL);
	_curveTree->setHeaderLabel("Curves");
	_curveTree->setSelectionMode(QAbstractItemView::ExtendedSelection);
	_curveTree->setContextMenuPolicy(Qt::CustomContextMenu);
	_curveTree->setMaximumSize(QSize(250, 500));
	*/

	_seqMenu = NULL;

	layout->addItem(treeout);

	_tabs = new QTabWidget(window);
	layout->addWidget(_tabs);
//	connect(_tabs, QTabWidget::currentChanged, this, &Main::tabChanged);
	
	setCentralWidget(window);

	_view = new StructureView(NULL);
	_view->setMain(this);
	_tabs->addTab(_view, "Structure");

	_seqView = new SequenceView(NULL, _fMaster);
	_seqView->setMain(this);
	_tabs->addTab(_seqView, "Sequence");

	_curveView = new CurveView(NULL);
	_fMaster->setCurveView(_curveView);
	_tabs->addTab(_curveView, "Graphs");

	_diff = new DiffDisplay(NULL, NULL);
	_diff->setMain(this);
//	_tabs->addTab(_diff, "Difference view");

	_couple = new StructureView(NULL);
	_coupleDisplay = new CoupleDisplay(NULL, _couple);

//	_tabs->addTab(_coupleDisplay, "Couple view");

	connect(_pdbTree, &QTreeWidget::customContextMenuRequested,
	        this, &Main::structureMenu);

	connect(_fMaster, &QTreeWidget::customContextMenuRequested,
	        this, &Main::fastaMenu);
	
	connect(_pdbTree, &QTreeWidget::itemSelectionChanged,
	        this, &Main::clickedStructure);
	
	connect(_diffTree, &QTreeWidget::currentItemChanged,
	        this, &Main::clickedDifference);

	makeMenu();
	show();
}

void Main::makeMenu()
{
	QMenu *file = menuBar()->addMenu(tr("&File"));
	QAction *act = file->addAction(tr("&Load structures"));
	connect(act, &QAction::triggered, this, &Main::loadStructures);
	act = file->addAction(tr("&Load fastas"));
	connect(act, &QAction::triggered, this, &Main::loadFastas);
	act = file->addAction(tr("&Load sequence metadata"));
	connect(act, &QAction::triggered, this, &Main::loadMetadata);
	act = file->addAction(tr("&Write protein fastas"));
	connect(act, &QAction::triggered, this, &Main::writeFastas);
	act = file->addAction(tr("&Write fasta subset"));
	connect(act, &QAction::triggered, this, &Main::writeSubset);
	act = file->addAction(tr("&Write calculated mutations"));
	connect(act, &QAction::triggered, this, &Main::writeMutations);
	act = file->addAction(tr("&Clear fastas"));
	connect(act, &QAction::triggered, _fMaster, &FastaMaster::clear);

	QMenu *structures = menuBar()->addMenu(tr("&Structure"));

	act = structures->addAction(tr("&Require mutation"));
	connect(act, &QAction::triggered, this, 
	        &Main::mutationWindow);
	act = structures->addAction(tr("&Highlight mutations"));
	connect(act, &QAction::triggered, _fMaster, 
	        &FastaMaster::highlightMutations);
	act = structures->addAction(tr("&Clear mutations"));
	connect(act, &QAction::triggered, _fMaster, 
	        &FastaMaster::clearMutations);
	act = structures->addAction(tr("&Generate sliding window"));
	connect(act, &QAction::triggered, this, &Main::prepareSlidingWindow);
	
	makeSequenceMenu();
}

void Main::mutationWindow()
{
	MutationWindow *w = new MutationWindow(NULL);
	w->setMain(this);
	w->show();
}

void Main::prepareSlidingWindow()
{
	_sw = new SlidingWindow(NULL, this);
	_sw->show();
}

void Main::slidingWindow(size_t window_size, std::string requirements, 
                         bool over)
{
	std::string folder = openDialogue(this, "Image folder", "", true, true);

	if (folder == "")
	{
		return;
	}
	
	if (_sw != NULL)
	{
		_sw->hide();
		_sw->deleteLater();
		_sw = NULL;
	}
	
	_fMaster->slidingWindowHighlight(_view, folder, window_size,
	                                 requirements, over);
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

void Main::loadFastas()
{
	LoadFastas *load = new LoadFastas(NULL);
	load->setMain(this);
	load->show();
}

void Main::loadMetadata()
{
	std::string filename = openDialogue(this, "Write sequence file", 
	                                    "Fasta file (*.fasta)", false);

	if (filename == "")
	{
		return;
	}
	
	_fMaster->loadMetadata(filename);
	makeSequenceMenu();
}

void Main::writeSubset()
{
	std::string filename = openDialogue(this, "Write sequence file", 
	                                    "Fasta file (*.fasta)", true);

	if (filename == "")
	{
		return;
	}
	
	_fMaster->writeOutFastas(filename);
}

void Main::writeMutations()
{
	std::string filename = openDialogue(this, "Write mutation metadata csv", 
	                                    "Fasta file (*.fasta)", true);
	
	_fMaster->writeOutMutations(filename);
}

void Main::writeFastas()
{
	std::string filename = openDialogue(this, "Write sequence file", 
	                                    "Fasta file (*.fasta)", true);

	if (filename == "")
	{
		return;
	}
	
	_fMaster->writeOutFastas(filename);
}

void Main::makeSequenceMenu()
{
	if (!_fMaster->isActive())
	{
		return;
	}

	if (_seqMenu == NULL)
	{
		_seqMenu = menuBar()->addMenu(tr("&Sequences"));
	}

	_fMaster->makeMenu(_seqMenu);
}

void Main::receiveSequence(Fasta *f)
{
	_ref->processNucleotides(f);
	_fMaster->addFasta(f);
}

void Main::receiveEnsemble(Ensemble *e)
{
	_pdbTree->addTopLevelItem(e);
	_view->addEnsemble(e);
	
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
	
	d->toCoupleDisplay(_coupleDisplay);
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
	
	for (size_t i = 0; i < e->chainCount(); i++)
	{
		std::cout << "Chain " << e->chain(i) <<
		": " << e->generateSequence(e->chain(i)) << std::endl;
	}

	_ref = e;
	_fMaster->setReference(e);
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

void Main::fastaMenu(const QPoint &p)
{
	QMenu *m = new QMenu();
	QPoint pos = centralWidget()->mapToGlobal(p);
	
	if (_fMaster->selectedItems().size() > 1)
	{
		_fMaster->makeGroupMenu(m);
		m->exec(pos);
		return;
	}

	FastaGroup *group = _fMaster->selectedGroup();

	if (group != NULL)
	{
		group->giveMenu(m);
	}

	if (_fMaster->selectedFasta())
	{
		m->addSeparator();
		_fMaster->selectedFasta()->giveMenu(m, group);
	}

	m->exec(pos);
}

void Main::tabChanged(int i)
{

}
