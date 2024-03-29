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

namespace QtCharts{
	class QChartView;
}
using namespace QtCharts;

class MyDictator;
class Database;
class DiffDisplay;
class CoupleDisplay;
class SlidingWindow;
class CurveView;
class QTabWidget;
class Ensemble;
class Fasta;
class FastaMaster;
class StructureView;
class SequenceView;
class QMenu;

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
	
	CoupleDisplay *coupleView()
	{
		return _coupleDisplay;
	}
	
	FastaMaster *fMaster()
	{
		return _fMaster;
	}

	void receiveEnsemble(Ensemble *e);
	void receiveSequence(Fasta *f);
	
	void makeReference(Ensemble *e);
	void setCommandLineArgs(int argc, char *argv[]);

	void makeSequenceMenu();
	
	void slidingWindow(size_t window_size, std::string requirements, 
	                   bool over);
public slots:
	void loadFastas();
	void writeFastas();
	void writeSubset();
	void loadMetadata();
	void scanMutation();
	void mutationWindow();
	void writeMutations();
	void loadStructures();
	void makeDifference();
	void updateDatabase();
	void getFromDatabase();
	void tabChanged(int i);
	void clickedStructure();
	void clickedDifference();
	void layeredScreenshot();
	void setChosenAsReference();
	void prepareSlidingWindow();
	void fastaMenu(const QPoint &p);
	void structureMenu(const QPoint &p);
protected:
	virtual void resizeEvent(QResizeEvent *e);
private:
	void makeMenu();
	QTreeWidget *_pdbTree;
	QTreeWidget *_diffTree;
	SlidingWindow *_sw;
	QTabWidget *_tabs;
	StructureView *_view;
	StructureView *_couple;
	SequenceView *_seqView;
	CoupleDisplay *_coupleDisplay;
	QChartView *_chartView;
	CurveView *_curveView;
	DiffDisplay *_diff;
	MyDictator *_dictator;
	QMenu *_seqMenu;
	
	Ensemble *_ref;
	FastaMaster *_fMaster;
	Database *_db;
	std::vector<std::string> _args;
};


#endif
