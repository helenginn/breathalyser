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

#include "LoadStructure.h"
#include "Ensemble.h"
#include "Main.h"

#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <hcsrc/FileReader.h>
#include <h3dsrc/Dialogue.h>
#include <libsrc/Multistate.h>
#include <libsrc/PDBReader.h>
#include <c4xsrc/QuickAtoms.h>

LoadStructure::LoadStructure(QWidget *parent) : QMainWindow(parent)
{
	setGeometry(500, 300, 540, 400);
	show();

	QLabel *l = new QLabel("PDB file:", this);
	l->setGeometry(20, 20, 100, 30);
	l->show();
	
	QLineEdit *e = new QLineEdit(this);
	e->setGeometry(120, 20, 300, 30);
	e->show();
	_pdbLine = e;
	
	QPushButton *b = new QPushButton("Choose...", this);
	b->setGeometry(420, 20, 100, 30);
	b->show();
	connect(b, &QPushButton::clicked, this, &LoadStructure::choosePDB);
	
	b = new QPushButton("Load", this);
	b->setGeometry(310, 350, 100, 30);
	b->show();

	connect(b, &QPushButton::clicked, this, &LoadStructure::loadPDB);
	
	b = new QPushButton("Done", this);
	b->setGeometry(420, 350, 100, 30);
	b->show();

	connect(b, &QPushButton::clicked, this, &LoadStructure::hide);
	connect(b, &QPushButton::clicked, this, &LoadStructure::deleteLater);
}

void LoadStructure::choosePDB()
{
	std::string filename = openDialogue(this, "Choose PDB", 
	                                    "Protein data bank (*.pdb)");

	_pdbLine->setText(QString::fromStdString(filename));
}

void LoadStructure::loadPDB()
{
	std::string pdb = _pdbLine->text().toStdString();
	std::string name = getBaseFilename(pdb);
	Multistate ms(pdb);
	ms.ignoreAtomsExcept("CA");
	ms.process();
	
	QuickAtoms *master = new QuickAtoms(NULL);
	Ensemble *top = new Ensemble(NULL, master);
	top->setName(name);
	
	if (ms.crystalCount() == 0)
	{
		PDBReader reader;
		reader.setFilename(pdb);
		reader.ignoreAtomsExcept("CA");
		CrystalPtr crystal = reader.getCrystal();
		
		if (crystal)
		{
			delete master;
			delete top;

			QuickAtoms *qa = new QuickAtoms(crystal);
			qa->fetchAtoms();
			master = qa;
			top = new Ensemble(NULL, master);
			top->setName(name);
		}
	}
	
	for (size_t i = 0; i < ms.crystalCount(); i++)
	{
		CrystalPtr crystal = ms.crystal(i);
		QuickAtoms *qa = new QuickAtoms(crystal);
		qa->fetchAtoms();
		Ensemble *conformer = new Ensemble(top, qa);
		master->addAtomsFrom(qa);
		conformer->setName(name + "_" + i_to_str(i));
	}

	if (ms.crystalCount() > 1)
	{
		master->divideThrough();
	}
	_main->receiveEnsemble(top);

	_pdbLine->setText("");
}
