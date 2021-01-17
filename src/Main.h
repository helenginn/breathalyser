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

#ifndef __breathalyser__main__
#define __breathalyser__main__

#include <QMainWindow>
#include <QTreeWidget>

class MyDictator;
class DiffDisplay;
class QTabWidget;
class Ensemble;
class StructureView;

class Main : public QMainWindow
{
Q_OBJECT
public:
	Main(QWidget *parent = NULL);
	
	int ensembleCount()
	{
		return _pdbTree->topLevelItemCount();
	}
	
	Ensemble *reference()
	{
		return _ref;
	}
	
	Ensemble *ensemble(int i);
	
	QTabWidget *tabs()
	{
		return _tabs;
	}
	
	StructureView *coupleView()
	{
		return _couple;
	}

	void receiveEnsemble(Ensemble *e);
	
	void makeReference(Ensemble *e);
	void setCommandLineArgs(int argc, char *argv[]);

public slots:
	void loadStructures();
	void setChosenAsReference();
	void clickedStructure();
	void clickedDifference();
	void structureMenu(const QPoint &p);
	void makeDifference();
protected:
	virtual void resizeEvent(QResizeEvent *e);
private:
	void makeMenu();
	QTreeWidget *_pdbTree;
	QTreeWidget *_diffTree;
	QTabWidget *_tabs;
	StructureView *_view;
	StructureView *_couple;
	DiffDisplay *_diff;
	MyDictator *_dictator;
	
	Ensemble *_ref;
	std::vector<std::string> _args;
};


#endif
