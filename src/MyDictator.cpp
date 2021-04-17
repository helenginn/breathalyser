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

#include "Main.h"
#include "MyDictator.h"
#include "FastaMaster.h"
#include "Fasta.h"
#include "LoadStructure.h"
#include "LoadFastas.h"
#include <iostream>
#include <hcsrc/FileReader.h>

MyDictator::MyDictator(Main *main) : Dictator()
{
	_main = main;
	_start = -1;
	_end = -1;
}

bool MyDictator::processRequest(std::string first, std::string last)
{
	if (first == "load-pdb")
	{
		std::vector<std::string> pdbs = split(last, ',');
		LoadStructure ls;
		ls.setMain(_main);
		
		for (size_t i = 0; i < pdbs.size(); i++)
		{
			ls.loadPDB(pdbs[i]);
		}
	}
	if (first == "focus-range")
	{
		std::vector<std::string> bits = split(last, '-');
		if (bits.size() == 2)
		{
			_start = atoi(bits[0].c_str());
			_end = atoi(bits[1].c_str());
		}
	}
	if (first == "load-nucleotide-seq")
	{
		std::vector<std::string> files = split(last, ',');
		LoadFastas lf;
		lf.setMain(_main);
		
		for (size_t i = 0; i < files.size(); i++)
		{
			lf.loadSequence(files[i], _start, _end, false);
		}
	}
	if (first == "load-protein-seq")
	{
		std::vector<std::string> files = split(last, ',');
		LoadFastas lf;
		lf.setMain(_main);
		
		for (size_t i = 0; i < files.size(); i++)
		{
			lf.loadSequence(files[i], _start, _end, true);
		}
	}
	if (first == "write-fastas")
	{
		_main->fMaster()->writeOutFastas(last);
	}
	if (first == "write-mutations")
	{
		_main->fMaster()->writeOutMutations(last);
	}
	if (first == "clear-fastas")
	{
		_main->fMaster()->clear();
	}
	if (first == "load-metadata")
	{
		_main->fMaster()->loadMetadata(last);
		_main->makeSequenceMenu();
	}
	if (first == "justify")
	{
		Fasta::setJustify(true);
	}
	if (first == "order-by")
	{
		_main->fMaster()->reorderBy(last);
	}
	if (first == "highlight-mutations")
	{
		_main->fMaster()->highlightMutations();
	}
	if (first == "clear-mutations")
	{
		_main->fMaster()->clearMutations();
	}
	if (first == "require-mutation")
	{
		_main->fMaster()->setTopAsCurrent();
		_main->fMaster()->requireMutation(last);
	}
	if (first == "update-database")
	{
		_main->updateDatabase();
	}
	if (first == "quit")
	{
		exit(0);
	}

	return true;
}
