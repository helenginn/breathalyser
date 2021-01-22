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
#include "FastaMaster.h"
#include "Ensemble.h"
#include "Fasta.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <iomanip>

#include <QMenu>

#include <hcsrc/FileReader.h>

FastaMaster::FastaMaster()
{
	_active = false;
	_req = INT_MAX;
	_aa = '\0';
}

void FastaMaster::writeCluster4xFile(std::string filename)
{
	std::ofstream csv;
	csv.open(filename);

	for (size_t i = 0; i < _subfastas.size() - 1; i++)
	{
		for (size_t j = i + 1; j < _subfastas.size(); j++)
		{
			_subfastas[i]->carefulCompareWithFasta(_subfastas[j]);
			double muts = _subfastas[i]->mutationCount();
			muts /= 4;

			double score = exp(-(muts * muts));
			if (score != score)
			{
				score = 0;
			}

			csv << _subfastas[i]->name() << ",";
			csv << _subfastas[j]->name() << ",";
			csv << std::setprecision(6) << std::fixed << std::showpoint;
			csv << score << std::endl;
		}

		_subfastas[i]->carefulCompareWithFasta(_fastas[0]);
	}

	csv.close();
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
		
		count++;

		muts << f->at(i)->name() << ",";
		muts << f->at(i)->mutationSummary() << std::endl;
	}
	
	muts.close();

	std::cout << "Written out mutations: " << count << " fastas of " 
	<< _fastas.size() << " total." << std::endl;
}

void FastaMaster::writeOutFastas(std::string filename, bool all)
{
	std::ofstream seqs;
	seqs.open(filename);
	
	int count = 0;
	
	std::vector<Fasta *> *f = &_fastas;
	
	if (!all)
	{
		f = &_subfastas;
	}

	for (size_t i = 0; i < f->size(); i++)
	{
		if (!f->at(i)->hasResult())
		{
			std::cout << "Skipping " << f->at(i)->name()
			<< " as could not locate protein sequence." << std::endl;
			continue;
		}
		
		count++;
		seqs << ">" << f->at(i)->name() << std::endl;
		seqs << f->at(i)->result() << std::endl;
	}
	
	seqs.close();

	std::cout << "Written out " << count << " fastas of " 
	<< _fastas.size() << " total." << std::endl;
}

void FastaMaster::addFasta(Fasta *f)
{
	_fastas.push_back(f);
	_names[f->name()] = f;
	
	if (_nameKeys.count(f->name()))
	{
		_keys[f] = _nameKeys[f->name()];
	}
	
	if (_fastas.size() > 1)
	{
		if (f->hasResult())
		{
			f->roughCompare(_fastas[0]->result(), 0);
		}

	}
	

	checkForMutation(f);

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
		Fasta *which = _names[components[0]];

		for (size_t j = 1; j < components.size(); j++)
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
		
		if (which == NULL)
		{
			continue;
		}
		
		_keys[which] = kv;

		count++;
	}
	
	std::cout << "Metadata for " << _nameKeys.size() << 
	" fastas in memory." << std::endl;
	std::cout << "Not assigned " << skip << " fastas not in memory." << std::endl;
	std::cout << "Loaded metadata for " << count << " fastas." << std::endl;
	std::cout << "Titles are: " << std::endl;
	
	_titles.reserve(_titles.size() + tmpTitles.size());
	_titles.insert(_titles.begin(), tmpTitles.begin(), tmpTitles.end());
	
	for (size_t i = 0; i < _titles.size(); i++)
	{
		std::cout << "\t" << _titles[i] << std::endl;
	}
	
	checkForMutations();
}

void FastaMaster::checkForMutation(Fasta *f)
{
	if (_keys[f].count("mutations"))
	{
		std::string val = _keys[f]["mutations"];
		f->loadMutations(val);
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

void FastaMaster::highlightRange(int start, int end)
{
	if (start >= end)
	{
		std::cout << "Start/end range invalid." << std::endl;
		return;
	}
	
	if (start < 0)
	{
		start = 0;
	}
	
	if (end > (int)_fastas.size())
	{
		end = _fastas.size();
	}

	_ref->clearBalls();
	_subfastas.clear();
	_subfastas.push_back(_fastas[0]);

	std::vector<std::string> mutations;
	int total = end - start;
	int count = 0;
	int step = total / 100;

	for (size_t j = start; j < (size_t)end; j++)
	{
		if (!_fastas[j]->hasResult())
		{
			continue;
		}

		if (_fastas[j] == _fastas[0])
		{
			continue;
		}

		if (!_fastas[j]->hasCompared())
		{
			_fastas[j]->carefulCompareWithFasta(_fastas[0]);

			if (count % step == 0)
			{
				std::cout << "." << std::flush;
			}
		}

		bool taken = _ref->processFasta(_fastas[j], _requirements);
		
		if (taken)
		{
			_subfastas.push_back(_fastas[j]);
		}
	}
	
	int seqCount = _ref->makeBalls();
	std::cout << "Number of sequences: " << seqCount << std::endl;
}

void FastaMaster::highlightMutations()
{
	std::cout << "Highlighting mutations." << std::endl;
	if (_fastas.size() == 0)
	{
		std::cout << "No fastas to highlight with." << std::endl;
	}
	std::cout << "Reference is " << _fastas[0]->name() << std::endl;
	highlightRange(0, _fastas.size());
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
		highlightRange(i, i + window);
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
			std::string value = _keys[_fastas[i]][_lastOrdered];
			key += value;
			
			int end = i + window;
			if (i + window > _fastas.size())
			{
				end = _fastas.size() - 1;
				value += " to end";
			}
			else
			{
				value += " to ";
				value += _keys[_fastas[end]][_lastOrdered];
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
	_fastas.clear();
	_lastOrdered = "";
}


std::string FastaMaster::valueForKey(Fasta *f, std::string key)
{
	if (_keys.count(f) == 0)
	{
		return "";
	}
	
	KeyValue kv = _keys[f];

	if (kv.count(key) == 0)
	{
		return "";
	}
	
	return kv[key];
}

typedef struct
{
	Fasta *f;
	std::string value;
} FastaValue;

bool smaller_value(const FastaValue &v1, const FastaValue &v2)
{
	return (v1.value < v2.value);
}

void FastaMaster::reorderBy(std::string title)
{
	if (std::find(_titles.begin(), _titles.end(), title) == _titles.end())
	{
		std::cout << "Cannot find title in database." << std::endl;
		return;
	}
	
	std::vector<FastaValue> values;

	for (size_t i = 0; i < _fastas.size(); i++)
	{
		std::string value = valueForKey(_fastas[i], title);
		FastaValue pair;
		pair.f = _fastas[i];
		pair.value = value;
		values.push_back(pair);
	}
	
	std::sort(values.begin(), values.end(), smaller_value);
	
	std::vector<Fasta *> better;

	for (size_t i = 0; i < values.size(); i++)
	{
		better.push_back(values[i].f);
	}
	
	_fastas = better;
	
	_lastOrdered = title;
	
	std::cout << "Reordered by " << title << "." << std::endl;
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
	_refSeq = e->generateSequence(_ref->chain(0), &_minRes);
}
