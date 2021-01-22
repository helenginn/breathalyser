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

#ifndef __breathalyser__fastamaster__
#define __breathalyser__fastamaster__

#include <vector>
#include <string>
#include <QTreeWidget>

class Fasta;
class QMenu;
class Ensemble;
class StructureView;

typedef std::map<std::string, std::string> KeyValue;
typedef std::map<Fasta *, KeyValue> FastaKeys;
typedef std::map<std::string, KeyValue> NameKeys;
typedef std::map<std::string, Fasta *> FastaNames;

class FastaMaster : public QTreeWidget
{
Q_OBJECT
public:
	FastaMaster(QWidget *parent);
	
	bool isActive()
	{
		return _active;
	}

	void addFasta(Fasta *f);
	
	size_t fastaCount()
	{
		return _fastas.size();
	}
	
	void setReference(Ensemble *e);

	void writeOutFastas(std::string filename, bool all = true);
	void writeOutMutations(std::string filename, bool all = true);
	void writeCluster4xFile(std::string filename);
	void loadMetadata(std::string metadata);
	void checkForMutations();
	void checkForMutation(Fasta *f);
	
	void reorderBy(std::string title);
	std::string valueForKey(Fasta *f, std::string key);
	void slidingWindowHighlight(StructureView *view,
	                            std::string folder, size_t window_size,
	                            std::string requirements, bool over);

	void makeMenu(QMenu *m);
	
	void requireMutation(std::string reqs)
	{
		_requirements = reqs;
	}
public slots:
	void highlightMutations();
	void clearMutations();
	void clear();
private:
	void highlightRange(int start, int end);

	int _req;
	int _minRes;
	unsigned char _aa;
	bool _active;
	Ensemble *_ref;
	std::string _refSeq;

	std::vector<std::string> _titles;
	std::vector<Fasta *> _fastas;
	std::vector<Fasta *> _subfastas;
	std::string _lastOrdered;
	std::string _requirements;
	
	FastaNames _names;
	FastaKeys _keys;
	NameKeys _nameKeys;
};

#endif
