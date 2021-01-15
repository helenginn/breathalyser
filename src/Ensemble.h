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

class QuickAtoms;

class Ensemble : public QTreeWidgetItem, public SlipObject
{
public:
	Ensemble(Ensemble *parent, QuickAtoms *qa);
	
	void setName(std::string name);

	QuickAtoms *getQuickAtoms()
	{
		return _qa;
	}

	void addCAlpha(vec3 point);
	void repopulate();
private:
	QuickAtoms *_qa;

	std::string _name;
};

#endif
