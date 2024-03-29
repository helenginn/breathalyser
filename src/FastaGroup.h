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

#ifndef __breathalyser__fastagroup__
#define __breathalyser__fastagroup__

#include <QTreeWidgetItem>
#include <vector>
#include <map>
#include <c4xsrc/Screen.h>

class QMenu;
class Fasta;
class Ensemble;
class FastaMaster;

class FastaGroup : public QObject, public QTreeWidgetItem, C4XAcceptor
{
Q_OBJECT
public:
	FastaGroup(FastaMaster *master);
	FastaGroup(FastaGroup *group);

	void initialise();
	
	void setEnsemble(Ensemble *e)
	{
		_ensemble = e;
	}

	void updateText();
	
	void setCustomName(std::string string)
	{
		_customName = string;
		updateText();
	}
	
	void setPermanent(bool p)
	{
		_permanent = p;
	}

	void setRequirements(std::string requirements)
	{
		_requirements = requirements;
	}

	void removeFasta(Fasta *f);
	void addFasta(Fasta *f);
	void addGroup(FastaGroup *g);
	
	size_t fastaCount()
	{
		return _fastas.size();
	}

	Fasta *fasta(int i)
	{
		return _fastas[i];
	}
	
	std::string lastOrdered()
	{
		return _lastOrdered;
	}
	
	void giveMenu(QMenu *m);
	void highlightOne(Fasta *f);
	void highlightRange(int start = 0, int end = 0);
	void writeOutFastas(std::string filename);
	void writeAlignments(std::string filename);

	typedef struct
	{
		Fasta *f;
		std::string value;
	} FastaValue;

	virtual void finished();

	void refreshToolTips();
	std::string generateText();
	std::string shortText();
	std::string countDescription();
	void fetchValues(std::string title);

	static void makeCurve(std::vector<FastaGroup *> groups,
	                      std::string title, std::string filename);
public slots:
	void split(std::string title, int bins, bool reorder);
	FastaGroup *makeRequirementGroup(std::string reqs);
	bool reorderBy(std::string title);
	void prepareCluster4x();
	void selectInverse();
	void removeGroup();
	void highlight();
protected:
	virtual void setData(int column, int role, const QVariant &value);
	void setModelData(QWidget *editor, QAbstractItemModel *model, 
                              const QModelIndex &index );

private:
	void addHighlight(Fasta *f);
	void countMutations();
	int lostMutations(size_t total);
	static bool smaller_value(const FastaGroup::FastaValue &v1, 
	                          const FastaGroup::FastaValue &v2);

	void titleLimits(double *min, double *max);
	int numberBetween(double min, double max);
	
	Screen *_screen;

	void clearFastas();

	std::vector<Fasta *> _fastas;
	std::string _requirements;

	Ensemble *_ensemble;
	FastaMaster *_master;
	FastaGroup *_group;
	
	std::map<std::string, Fasta *> _nameMap;

	bool _permanent;
	std::string _lastOrdered;

	std::string _customName;
	std::map<std::string, int> _mutCounts;
	std::vector<std::string> _muts;
};

#endif
