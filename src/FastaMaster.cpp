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

#include "StructureView.h"
#include "SequenceView.h"
#include "WidgetFasta.h"
#include "FastaMaster.h"
#include "FastaGroup.h"
#include "Ensemble.h"
#include "Fasta.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <iomanip>

#include <QMenu>

#include <hcsrc/FileReader.h>
#include <h3dsrc/Dialogue.h>

FastaMaster::FastaMaster(QWidget *parent) : QTreeWidget(parent)
{
	_seqView = NULL;
	_top = new FastaGroup(this);
	_top->setPermanent(true);
	addTopLevelItem(_top);
	_top->setCustomName("all");
	setCurrentItem(_top);

	_active = false;
	_req = INT_MAX;
	_aa = '\0';
	
	connect(this, &QTreeWidget::itemClicked, this,
	        &FastaMaster::itemClicked);
}

void FastaMaster::writeOutMutations(std::string filename, bool all)
{
	if (filename.length() == 0)
	{
		std::cout << "No filename." << std::endl;
		exit(0);
	}

	std::ofstream muts;
	muts.open(filename);
	
	int count = 0;
	
	std::vector<Fasta *> *f = &_fastas;
	
	if (!all)
	{
		f = &_subfastas;
	}
	
	muts << "sequence_name,mutations" << std::endl;

	for (size_t i = 0; i < f->size(); i++)
	{
		if (!f->at(i)->hasCompared())
		{
			f->at(i)->carefulCompareWithFasta(_fastas[0]);
		}
		
		if (f->at(i)->isProblematic())
		{
			continue;
		}
		
		count++;

		muts << f->at(i)->name() << ",";
		muts << f->at(i)->mutationSummary() << std::endl;
	}
	
	muts.close();

	std::cout << "Written out mutations: " << count << " fastas of " 
	<< _fastas.size() << " total." << std::endl;
}

void FastaMaster::writeOutFastas(std::string filename)
{
	if (selectedGroup())
	{
		selectedGroup()->writeOutFastas(filename);
	}
}

void FastaMaster::addFasta(Fasta *f)
{
	_fastas.push_back(f);
	_names[f->name()] = f;
	
	if (_fastas.size() > 1)
	{
		if (f->hasResult())
		{
			f->roughCompare(_fastas[0]->result(), 0);
		}

	}
	
	checkForMutation(f);
	_top->addFasta(f);
	_top->updateText();

	_active = true;
}

void FastaMaster::loadMetadata(std::string fMetadata)
{
	if (!file_exists(fMetadata))
	{
		std::cout << "File does not exist" << std::endl;
		return;
	}

	std::string contents = get_file_contents(fMetadata);
	std::vector<std::string> lines = split(contents, '\n');
	
	if (lines.size() == 0)
	{
		std::cout << "File empty" << std::endl;
		return;
	}
	
	std::string header = lines[0];
	std::vector<std::string> titles = split(header, ',');

	std::vector<std::string> tmpTitles;
	
	for (size_t i = 0; i < titles.size(); i++)
	{
		tmpTitles.push_back(titles[i]);
	}
	
	if (tmpTitles.size() == 0)
	{
		std::cout << "Lines are empty" << std::endl;
		return;
	}
	
	std::cout << "We assume " << tmpTitles[0] << " is the "\
	"sequence identifier. If this is not the case, "\
	"fix and reload" << std::endl;
	
	size_t skip = 0;
	size_t count = 0;
	
	for (size_t i = 1; i < lines.size(); i++)
	{
		std::vector<std::string> components = split(lines[i], ',');
		
		if (components.size() != tmpTitles.size())
		{
			std::cout << "Skipping line " << i << " - incorrect number "\
			" of entries." << std::endl;
			std::cout << "\t" << lines[i] << std::endl;
		}
		
		KeyValue kv;

		trim(components[0]);
		
		if (_nameKeys.count(components[0]))
		{
			kv = _nameKeys[components[0]];
		}

		for (size_t j = 0; j < components.size(); j++)
		{
			trim(components[j]);
			kv[tmpTitles[j]] = components[j];
		}

		_nameKeys[components[0]] = kv;
		
		if (_names.count(components[0]) == 0)
		{
			skip++;
			continue;
		}

		count++;
	}
	
	std::cout << "Metadata for " << _nameKeys.size() << 
	" fastas in memory." << std::endl;
	std::cout << "Not assigned " << skip << " fastas not in memory." << std::endl;
	std::cout << "Loaded metadata for " << count << " fastas." << std::endl;
	std::cout << "Titles are: " << std::endl;
	
	_titles.reserve(_titles.size() + tmpTitles.size());
	
	for (size_t i = 0; i < tmpTitles.size(); i++)
	{
		if (!hasKey(tmpTitles[i]))
		{
			_titles.push_back(tmpTitles[i]);
		}
	}
	
	for (size_t i = 0; i < _titles.size(); i++)
	{
		std::cout << "\t" << _titles[i] << std::endl;
	}
	
	checkForMutations();
}

void FastaMaster::requireMutation(std::string reqs)
{
	if (selectedGroup() == NULL)
	{
		setTopAsCurrent();
	}

	selectedGroup()->makeRequirementGroup(reqs);
}

void FastaMaster::checkForMutation(Fasta *f)
{
	if (fastaCount() == 0)
	{
		return;
	}

	if (_nameKeys[f->name()].count("mutations"))
	{
		std::string val = _nameKeys[f->name()]["mutations"];
		f->loadMutations(val, fasta(0)->result());
		
		std::string new_val = f->mutationSummary();
		_nameKeys[f->name()]["mutations"] = new_val;
	}
}

void FastaMaster::checkForMutations()
{
	if (std::find(_titles.begin(), _titles.end(), "mutations") 
	    != _titles.end())
	{
		for (size_t i = 0; i < _fastas.size(); i++)
		{
			checkForMutation(_fastas[i]);
		}
	}
	
	std::cout << "Loaded mutations from metadata." << std::endl;
}

void FastaMaster::clearMutations()
{
	_ref->clearBalls();
}


void FastaMaster::highlightMutations()
{
	std::cout << "Highlighting mutations." << std::endl;

	if (_fastas.size() == 0)
	{
		std::cout << "No fastas to highlight with." << std::endl;
	}
	
	if (selectedGroup())
	{
		selectedGroup()->highlightRange();
	}
}

void FastaMaster::slidingWindowHighlight(StructureView *view,
                                         std::string folder, size_t window,
                                         std::string requirements, bool over)
{
	if (requirements.length())
	{
		_requirements = requirements;
	}

	size_t step = window / 10;
	if (step == 0)
	{
		step = 1;
	}
	int count = 0;
	
	if (over)
	{
		std::string pattern = folder + "/fr*.png";
		std::vector<std::string> files = glob(pattern);
		
		for (size_t i = 0; i < files.size(); i++)
		{
			remove(files[i].c_str());
		}
	}
	
	FileReader::setOutputDirectory(folder);

	for (size_t i = 0; i < _fastas.size(); i += step)
	{
		if (_fastas.size() - i < window / 2)
		{
			continue;
		}

		std::cout << "Highlighting range " << std::endl;
		_top->highlightRange(i, i + window);
		view->update();
		
		std::string number = i_to_str(count);
		count++;
		std::string zeros;

		if (count <= 99999)
		{
			zeros = std::string(5 - number.length(), '0');
		}
		
		std::string key = "";
		
		if (_lastOrdered.length())
		{
			key = _lastOrdered + "_";
			std::string value = valueForKey(_fastas[i], _lastOrdered);
			key += value;
			
			if (i + window > _fastas.size())
			{
				value += " to end";
			}
			else
			{
				value += " to ";
				value += valueForKey(_fastas[i], _lastOrdered);
			}


			view->addLabel(value);
		}
		
		std::string filename = "fr_" + key + "_" + zeros + number + ".png";
		std::string path = FileReader::addOutputDirectory(filename);
		std::cout << path << std::endl;

		view->saveImage(path);
	}
	
	_req = INT_MAX;
	_aa = '\0';
}

void FastaMaster::clear()
{
	clear();

	for (size_t i = 0; i < _fastas.size(); i++)
	{
		_top->removeChild(_fastas[i]);
		delete _fastas[i];
	}

	addTopLevelItem(_top);
	_top->updateText();

	_fastas.clear();
	_lastOrdered = "";
}


void FastaMaster::reorderBy(std::string title)
{
	_top->reorderBy(title);
}

size_t FastaMaster::fastaCount()
{
	return _top->fastaCount();
}

Fasta *FastaMaster::fasta(int i)
{
	return _top->fasta(i);
}

void FastaMaster::makeGroupMenu(QMenu *m)
{
	std::cout << "Bleh" << std::endl;
	QMenu *submenu = m->addMenu(tr("&Relative population CSV..."));

	for (size_t i = 0; i < _titles.size(); i++)
	{
		QString qTitle = QString::fromStdString(_titles[i]);

		QAction *act = submenu->addAction(qTitle);
		connect(act, &QAction::triggered, 
		        this, [=]() { makeCurves(_titles[i]); });
	}
}

void FastaMaster::makeMenu(QMenu *m)
{
	m->clear();

	QMenu *submenu = m->addMenu(tr("&Reorder by..."));
	
	for (size_t i = 0; i < _titles.size(); i++)
	{
		QString qTitle = QString::fromStdString(_titles[i]);

		QAction *act = submenu->addAction(qTitle);
		connect(act, &QAction::triggered, 
		        this, [=]() {reorderBy(_titles[i]);});
	}

	std::cout << "Currently loaded: " << _fastas.size() << " fastas." << std::endl;
}

void FastaMaster::setReference(Ensemble *e)
{
	_ref = e;
	_top->setEnsemble(_ref);
	_refSeq = e->generateSequence(_ref->chain(0), &_minRes);
}

std::vector<FastaGroup *> FastaMaster::selectedGroups()
{
	std::vector<FastaGroup *> groups;
	QList<QTreeWidgetItem *> items = selectedItems();
	
	for (int i = 0; i < items.size(); i++)
	{
		QTreeWidgetItem *item = items[i];
		FastaGroup *g = dynamic_cast<FastaGroup *>(item);
		if (g != NULL)
		{
			groups.push_back(g);
		}
	}

	return groups;
}

FastaGroup *FastaMaster::selectedGroup()
{
	QTreeWidgetItem *item = currentItem();
	WidgetFasta *f = dynamic_cast<WidgetFasta *>(item);
	
	if (f != NULL)
	{
		return f->group();
	}

	FastaGroup *g = dynamic_cast<FastaGroup *>(item);

	return g;
}

Fasta *FastaMaster::selectedFasta()
{
	QTreeWidgetItem *item = currentItem();
	WidgetFasta *f = dynamic_cast<WidgetFasta *>(item);
	
	if (f == NULL)
	{
		return NULL;
	}

	return f->fasta();
}

void FastaMaster::itemClicked(QTreeWidgetItem *item, int column)
{
	FastaGroup *grp = selectedGroup();
	if (grp != NULL)
	{
		grp->highlightRange();
		_seqView->populate(grp);
	}

	Fasta *f = selectedFasta();
	
	if (f != NULL && _seqView != NULL)
	{
		_seqView->populate(f);
	}
}

bool FastaMaster::fastaHasKey(Fasta *f, std::string key)
{
	if (_nameKeys.count(f->name()) == 0)
	{
		return false;
	}
	
	if (_nameKeys[f->name()].count(key) == 0)
	{
		return false;
	}
	
	return true;
}

std::string FastaMaster::valueForKey(Fasta *f, std::string key)
{
	if (fastaHasKey(f, key))
	{
		return _nameKeys[f->name()][key];
	}
	
	return "";
}

bool FastaMaster::hasKey(std::string key)
{
	return (std::find(_titles.begin(), _titles.end(), key) != _titles.end());
}

bool FastaMaster::isReference(Fasta *f)
{
	return _top->fastaCount() && f == _top->fasta(0);
}

void FastaMaster::makeCurves(std::string title)
{
	std::string filename = openDialogue(this, "Write CSV file", 
	                                    "Comma separated values (*.csv)", 
	                                    false);

	makeCurvesForFilename(filename, title);
}

void FastaMaster::makeCurvesForFilename(std::string filename, 
                                        std::string title)
{
	std::vector<FastaGroup *> groups = selectedGroups();
	std::cout << groups.size() << std::endl;
	FastaGroup::makeCurve(groups, title, filename);
}

void FastaMaster::setTopAsCurrent()
{
	setCurrentItem(_top);
}

