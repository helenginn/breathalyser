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

#ifndef __breathalyser__fasta__
#define __breathalyser__fasta__

#include <string>
#include <map>
#include <vector>

#include <QTreeWidgetItem>

class FastaGroup;

typedef std::map<int, int> IntMap;

class Fasta : public QObject, public QTreeWidgetItem
{
Q_OBJECT
public:
	Fasta(std::string name);
	
	void setSequence(std::string seq, bool protein);
	void description();

	bool findNextORF();
	bool findNextATG();
	bool findNextStopCodon();
	
	void setOffset(int off)
	{
		_offset = off;
	}

	std::string name()
	{
		return _name;
	}
	
	std::string result()
	{
		return _result;
	}
	
	unsigned char letter(int i);
	
	bool hasResult()
	{
		return _result.length() > 0;
	}

	int gapCount();
	
	std::string mutation(int i)
	{
		return _mutations[i];
	}
	
	size_t mutationCount()
	{
		return _mutations.size();
	}
	
	bool hasMutation(std::string mut);
	
	bool hasCompared()
	{
		return _compared;
	}
	
	bool isProblematic()
	{
		return _problematic;
	}

	void carefulCompareWithFasta(Fasta *f);
	void carefulCompareWithString(std::string seq2);
	double compareWithFasta(Fasta *f);
	void addMutation(char fromWhat, int mut, char towhat);

	int oneSidedMutations(Fasta *other);
	int sharedMutations(Fasta *other);

	std::string mutationSummary();
	void sortMutations();

	std::string generateSequence();
	bool roughlyAlign(std::string mine, std::string ref, int minRes);
	std::string roughCompare(std::string seq, int minRes);
	void loadMutations(std::string muts, std::string ref);

	void leftJustifyDeletions();

	void giveMenu(QMenu *m, FastaGroup *g);
	
	void setLastValue(std::string v)
	{
		_lastValue = v;
	}
	
	std::string lastValue()
	{
		return _lastValue;
	}
signals:
	void refreshMutations();
public slots:
	void setIsProblematic();
private:
	void refreshToolTips();
	void organiseMap();
	std::string deletionSequence(int start, int end, int go_back);
	void nudgeMap(int start, int dir);

	bool _compared;
	bool _problematic;
	int _orf;
	int _offset;
	int _stop;
	std::string _name;
	std::string _seq;
	std::string _result;
	std::string _ref;
	std::string _lastValue;
	
	std::vector<std::string> _mutations;
	IntMap _meToRef;
	IntMap _refToMe;
};

#endif
