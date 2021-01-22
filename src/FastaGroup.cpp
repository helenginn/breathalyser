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

#include <algorithm>

#include "FastaGroup.h"
#include "FastaMaster.h"
#include "WidgetFasta.h"
#include "Fasta.h"
#include "Ensemble.h"
#include <QMenu>
#include <QStyledItemDelegate>
#include <hcsrc/FileReader.h>
#include <c4xsrc/ClusterList.h>
#include <c4xsrc/AveCSV.h>
#include <c4xsrc/Group.h>
#include <c4xsrc/Screen.h>

void FastaGroup::initialise()
{
	_permanent = false;
	_ensemble = _master->getReference();

	Qt::ItemFlags fl = flags();
	setFlags(fl | Qt::ItemIsEditable);
}

FastaGroup::FastaGroup(FastaMaster *master) : QTreeWidgetItem(master)
{
	_screen = NULL;
	_group = NULL;
	_master = master;
	initialise();
}

FastaGroup::FastaGroup(FastaGroup *group) : QTreeWidgetItem(group)
{
	_screen = NULL;
	_group = group;
	_master = group->_master;
	initialise();
}

std::string FastaGroup::generateText()
{
	std::string text;
	if (_customName.length() == 0)
	{
		text = "Group ";

		if (_requirements.length())
		{
			text += "requiring " + _requirements + " ";
		}
	}
	else
	{
		text = _customName + " ";
	}
	
	text += "(" + i_to_str(childCount()) + ")";
	return text;
}

void FastaGroup::setData(int column, int role, const QVariant &value)
{
	if (role == Qt::EditRole)
	{
		_customName = value.toString().toStdString();
	}
	
	QTreeWidgetItem::setData(column, role, value);
	emitDataChanged();
}

void FastaGroup::updateText()
{
	std::string text = generateText();
	setText(0, QString::fromStdString(text));
}

void FastaGroup::addFasta(Fasta *f)
{
	_fastas.push_back(f);
	WidgetFasta *item = new WidgetFasta(f, this);
	addChild(item);
	_nameMap[f->name()] = f;
	connect(f, &Fasta::refreshMutations, 
	        this, &FastaGroup::highlight);
}

void FastaGroup::removeFasta(Fasta *f)
{
	for (size_t i = 0; i < _fastas.size(); i++)
	{
		if (_fastas[i] == f)
		{
			disconnect(f, &Fasta::refreshMutations, 
			           this, &FastaGroup::highlight);

			_fastas.erase(_fastas.begin() + i);
		}
	}
}

void FastaGroup::addGroup(FastaGroup *g)
{
	if (g->QTreeWidgetItem::parent() == this)
	{
		removeChild(g);
	}

	insertChild(0, g);
}

void FastaGroup::highlight()
{
	highlightRange(0, 0);
}

void FastaGroup::highlightRange(int start, int end)
{
	if (start < 0)
	{
		start = 0;
	}
	
	if (end == 0 || end > (int)fastaCount())
	{
		end = fastaCount();
	}

	if (start > end)
	{
		std::cout << "Start/end range invalid." << std::endl;
		return;
	}
	
	_ensemble->clearBalls();

	int total = end - start;
	int max_stages = 100;
	int per_stage = total / max_stages;
	int stages = 0;
	int count = 0;
	
	if (fastaCount() == 0)
	{
		return;
	}

	std::cout << "Reference is " << fasta(0)->name() << std::endl;
	for (size_t j = start; j < (size_t)end; j++)
	{
		if (fasta(j) == fasta(0))
		{
			continue;
		}
		
		if (!fasta(j)->hasCompared())
		{
			fasta(j)->carefulCompareWithFasta(fasta(0));

			count++;
		}

		if (count > stages + per_stage)
		{
			std::cout << "." << std::flush;
		}

		_ensemble->processFasta(fasta(j), _requirements);
	}
	
	updateText();
	
	int seqCount = _ensemble->makeBalls();
	std::cout << "No. sequences displayed: " << seqCount << std::endl;
}

void FastaGroup::makeRequirementGroup(std::string reqs)
{
	if (fastaCount() <= 1)
	{
		return;
	}

	FastaGroup *g = _master->selectedGroup();
	FastaGroup *grp = NULL;
	
	bool isTop = (g == _master->topGroup());
	if (isTop)
	{
		grp = new FastaGroup(_master);
	}
	else
	{
		grp = new FastaGroup(this);
	}

	Fasta *ref = fasta(0);
	grp->setRequirements(reqs);
	grp->addFasta(ref);

	for (size_t i = 1; i < fastaCount(); i++)
	{
		if (_ensemble->shouldProcess(fasta(i), reqs))
		{
			grp->addFasta(fasta(i));
		}
	}
	
	grp->updateText();

	if (isTop)
	{
		_master->addTopLevelItem(grp);
	}
	else
	{
		addGroup(grp);
	}

	_master->setCurrentItem(grp);
	grp->highlightRange();
}

void FastaGroup::clearFastas()
{
	for (int i = 0; i < childCount(); i++)
	{
		if (dynamic_cast<WidgetFasta *>(child(i)) != NULL)
		{
			takeChild(i);
			i--;
		}
	}

	for (size_t i = 0; i < fastaCount(); i++)
	{
		removeFasta(fasta(i));
		i--;
	}
}

bool FastaGroup::smaller_value(const FastaGroup::FastaValue &v1, 
                               const FastaGroup::FastaValue &v2)
{
	return (v1.value < v2.value);
}

void FastaGroup::selectInverse()
{
	Fasta *ref = fasta(0);
	std::string name = "Not " + generateText();
	FastaGroup *grp = new FastaGroup(_master);
	grp->addFasta(ref);
	grp->setCustomName(name);

	for (size_t i = 1; i < _master->fastaCount(); i++)
	{
		Fasta *trial = _master->fasta(i);
		bool found = (std::find(_fastas.begin(), _fastas.end(), trial)
		              != _fastas.end());
		
		if (!found)
		{
			grp->addFasta(trial);
		}
	}

	_master->addTopLevelItem(grp);
}

void FastaGroup::giveMenu(QMenu *m)
{
	if (!_permanent)
	{
		QAction *act = m->addAction("Select inverse");
		connect(act, &QAction::triggered, this, &FastaGroup::selectInverse);
	}

	QMenu *submenu = m->addMenu(tr("&Reorder by..."));

	for (size_t i = 0; i < _master->titleCount(); i++)
	{
		std::string title = _master->title(i);
		QString qTitle = QString::fromStdString(title);

		QAction *act = submenu->addAction(qTitle);
		connect(act, &QAction::triggered, this, [=]() {reorderBy(title);});
	}
	
	QAction *act = m->addAction("Split by reordered");
	connect(act, &QAction::triggered, 
	        this, [=] { split("", 0, true); });
	
	act = m->addAction("Split using cluster4x");
	connect(act, &QAction::triggered, this, &FastaGroup::prepareCluster4x);

	if (!_permanent)
	{
		QAction *act = m->addAction("Remove group");
		connect(act, &QAction::triggered, this, &FastaGroup::removeGroup);
	}
}

void FastaGroup::removeGroup()
{
	if (_group == NULL)
	{
		int index = treeWidget()->indexOfTopLevelItem(this);
		treeWidget()->takeTopLevelItem(index);
	}
	else if (_group != NULL)
	{
		_group->removeChild(this);
	}
}

void FastaGroup::prepareCluster4x()
{
	AveCSV *csv = new AveCSV(NULL, "");
	csv->startNewCSV("Sequence similarity");
	
	std::vector<Fasta *> copy = _fastas;
	
	std::random_shuffle(copy.begin(), copy.end());
	
	int limit = 1000 + 1;
	if (limit > copy.size())
	{
		limit = copy.size();
	}

	for (size_t i = 1; i < limit - 1; i++)
	{
		for (size_t j = i + 1; j < limit; j++)
		{
			double muts = copy[i]->sharedMutations(copy[j]);
			muts /= 4;

			double score = exp(-(muts * muts));
			if (score != score)
			{
				score = 0;
			}

			csv->addValue(copy[i]->name(), copy[j]->name(), score);
			csv->addValue(copy[j]->name(), copy[i]->name(), score);
		}
	}

	if (_screen != NULL)
	{
		_screen->hide();
		_screen->deleteLater();
	}
	_screen = new Screen(NULL);
	_screen->setWindowTitle("cluster4x - sequins");
	_screen->setReturnJourney(this);
	ClusterList *list = _screen->getList();
	csv->setList(list);
	csv->preparePaths();
	csv->setChosen(0);

	Group *top = Group::topGroup();
	top->setCustomName(generateText());
	top->updateText();

	_screen->show();
}

void FastaGroup::finished()
{
	ClusterList *list = _screen->getList();

	for (size_t i = 0; i < list->groupCount(); i++)
	{
		Group *g = list->group(i);
		
		if (!g->isMarked())
		{
			continue;
		}

		FastaGroup *grp = new FastaGroup(this);
		Fasta *ref = fasta(0);
		grp->addFasta(ref);
		std::string def = "Split: " + generateText();
		
		std::string custom = g->getCustomName();
		if (custom.length())
		{
			def = custom;
		}

		for (size_t j = 0; j < g->mtzCount(); j++)
		{
			std::string fastaName = g->getMetadata(j);
			Fasta *myF = _nameMap[fastaName];
			grp->addFasta(myF);
		}

		grp->setCustomName(def);
		
		addGroup(grp);
		_master->setCurrentItem(grp);
	}
	
	_screen->hide();
}

void FastaGroup::refreshToolTips()
{
	for (size_t i = 0; i < childCount(); i++)
	{
		WidgetFasta *wf = dynamic_cast<WidgetFasta *>(child(i));
		
		if (wf)
		{
			wf->refreshTip();
		}
	}
}

void FastaGroup::writeOutFastas(std::string filename)
{
	std::ofstream seqs;
	seqs.open(filename);
	
	int count = 0;
	
	for (size_t i = 0; i < fastaCount(); i++)
	{
		if (!fasta(i)->hasResult())
		{
			std::cout << "Skipping " << fasta(i)->name()
			<< " as could not locate protein sequence." << std::endl;
			continue;
		}
		
		count++;
		seqs << ">" << fasta(i)->name() << std::endl;
		seqs << fasta(i)->result() << std::endl;
	}
	
	seqs.close();

	std::cout << "Written out " << count << " fastas of " 
	<< fastaCount() << " total." << std::endl;
}

void FastaGroup::reorderBy(std::string title)
{
	if (!_master->hasKey(title))
	{
		std::cout << "Cannot find title (" << title << 
		") in database." << std::endl;
		return;
	}
	
	if (fastaCount() == 0)
	{
		std::cout << "Nothing in group." << std::endl;
		return;
	}

	Fasta *ref = fasta(0);
	std::vector<FastaGroup::FastaValue> values;

	for (size_t i = 1; i < fastaCount(); i++)
	{
		std::string value = _master->valueForKey(fasta(i), title);
		FastaGroup::FastaValue pair;
		pair.f = fasta(i);
		pair.value = value;
		values.push_back(pair);
		
		fasta(i)->setLastValue(value);
	}
	
	std::sort(values.begin(), values.end(), FastaGroup::smaller_value);
	
	clearFastas();

	addFasta(ref);

	for (size_t i = 0; i < values.size(); i++)
	{
		addFasta(values[i].f);
	}
	
	_lastOrdered = title;
	refreshToolTips();
}

void FastaGroup::split(std::string title, int bins, bool reorder)
{
	if (title == "")
	{
		title = _lastOrdered;
	}
	
	if (bins == 0)
	{
		bins = fastaCount() / 1000 + 1;
	}

	if (!_master->hasKey(title))
	{
		std::cout << "Cannot find title in database." << std::endl;
		return;
	}

	if (reorder)
	{
		reorderBy(title);
	}
	
	int in_each = fastaCount() / bins + 1;
	
	for (size_t i = 0; i < fastaCount(); i += in_each)
	{
		FastaGroup *grp = new FastaGroup(this);
		Fasta *ref = fasta(0);
		grp->addFasta(ref);
		
		for (size_t j = i; j < i + in_each && j < fastaCount(); j++)
		{
			if (fasta(j) == ref)
			{
				continue;
			}

			grp->addFasta(fasta(j));
		}

		std::string def = "Binned: " + generateText();
		grp->setCustomName(def);
		
		addGroup(grp);
		_master->setCurrentItem(grp);
	}
}

