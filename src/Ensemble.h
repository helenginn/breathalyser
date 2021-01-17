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

#ifndef __breathalyser__ensemble__ 
#define __breathalyser__ensemble__ 

#include <QTreeWidgetItem>
#include <h3dsrc/SlipObject.h>

/******************** Crystal definition ***************************/
namespace Vagabond
{
	class Crystal;
}

#define ToCrystalPtr(a) (boost::static_pointer_cast<Vagabond::Crystal>((a)))
typedef boost::shared_ptr<Vagabond::Crystal> CrystalPtr;
typedef boost::weak_ptr<Vagabond::Crystal> CrystalWkr;

/****************** end Crystal definition *************************/

class Ensemble : public QTreeWidgetItem, public SlipObject
{
public:
	Ensemble(Ensemble *parent, CrystalPtr c);
	
	void setName(std::string name);
	
	std::string name()
	{
		return _name;
	}
	
	size_t chainCount()
	{
		return _chains.size();
	}
	
	std::string chain(int i)
	{
		return _chains[i];
	}

	CrystalPtr crystal()
	{
		return _crystal;
	}
	
	void setReference(bool ref)
	{
		_isReference = ref;
		updateText();
	}

	void addCAlpha(vec3 point);
	void repopulate();
	void updateText();
	vec3 averagePos();

	virtual void render(SlipGL *gl);

	void minMaxResidues(std::string ch, int *min, int *max);
	std::string generateSequence(std::string chain);
	std::string findMatchingChain(std::string ch, Ensemble *other);
private:
	void findChains();
	void convertToBezier();
	CrystalPtr _crystal;

	bool _isReference;
	std::string _name;
	std::vector<std::string> _chains;
};

#endif
