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

#include "SequenceView.h"
#include "FastaMaster.h"
#include "FastaGroup.h"
#include "Fasta.h"

#include <QScrollArea>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QCheckBox>
#include <QVBoxLayout>
#include <hcsrc/FileReader.h>

QBoxLayout *SequenceView::setupWindow()
{
	QWidget *w = centralWidget();
	if (w != NULL)
	{
		w->hide();
		w->deleteLater();
	}

	QScrollArea *area = new QScrollArea(NULL);
	QVBoxLayout *areaBox = new QVBoxLayout();
	area->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOn);
	area->setLayout(areaBox);
	setCentralWidget(area);
	
	return areaBox;
}

SequenceView::SequenceView(QWidget *parent, FastaMaster *m) 
: QMainWindow(parent)
{
	setStyleSheet("background-color:light gray;");
	setAutoFillBackground( true );

	_master = m;
	_master->setSequenceView(this);
	setupWindow();
	_fasta = NULL;
	_group = NULL;
}

void standardiseWidths(std::vector<QLabel *> labels)
{
	int biggest = 0;
	
	for (size_t i = 0; i < labels.size(); i++)
	{
		QLabel *l = labels[i];
		int width = l->fontMetrics().boundingRect(l->text()).width();
		biggest = std::max(width, biggest);
	}
	
	for (size_t i = 0; i < labels.size(); i++)
	{
		labels[i]->setMaximumSize(biggest + 10, 20);
	}
}

QLabel *SequenceView::addToLayout(QBoxLayout *box, std::string left, 
                                  std::string right)
{
	QHBoxLayout *hbox = new QHBoxLayout();
	hbox->setAlignment(Qt::AlignTop);
	hbox->setSizeConstraint(QLayout::SetMinAndMaxSize);
	QLabel *l = new QLabel(QString::fromStdString(left), this);
	QLabel *r = new QLabel(QString::fromStdString(right), this);
	hbox->addWidget(l);
	hbox->addWidget(r);
	box->addLayout(hbox);

	return l;
}

void SequenceView::populate(FastaGroup *g)
{
	_group = g;
	_group = NULL;
	QBoxLayout *box = setupWindow();

	std::vector<QLabel *> labels;
	QVBoxLayout *metaBox = new QVBoxLayout();
	metaBox->setAlignment(Qt::AlignTop);

	{
		std::string text = g->generateText();
		QLabel *l = addToLayout(metaBox, "Group name:", text);
		labels.push_back(l);
	}

	{
		std::string text = i_to_str(g->fastaCount());
		QLabel *l = addToLayout(metaBox, "Sequence count:", text);
		labels.push_back(l);
	}

	std::string text = g->countDescription();
	QLabel *l = addToLayout(metaBox, "Representative mutations:", text);
	labels.push_back(l);
	
	standardiseWidths(labels);
	box->addLayout(metaBox);
}

void SequenceView::populate(Fasta *f)
{
	_fasta = f;
	_group = NULL;
	QBoxLayout *box = setupWindow();

	std::vector<QLabel *> labels;
	QVBoxLayout *metaBox = new QVBoxLayout();
	metaBox->setAlignment(Qt::AlignTop);
	metaBox->setSizeConstraint(QLayout::SetMinAndMaxSize);

	for (size_t i = 0; i < _master->titleCount(); i++)
	{
		QHBoxLayout *hbox = new QHBoxLayout();
		QString title = QString::fromStdString(_master->title(i));
		std::string val = _master->valueForKey(f, _master->title(i));
		
		QLabel *l = new QLabel(title, this);
		labels.push_back(l);
		
		QLabel *r = new QLabel(QString::fromStdString(val), this);

		hbox->addWidget(l);
		hbox->addWidget(r);

		metaBox->addLayout(hbox);
	}
	
	standardiseWidths(labels);
	box->addLayout(metaBox);
}
