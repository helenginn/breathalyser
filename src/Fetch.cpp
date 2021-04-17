// splitseq
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

#include "FastaMaster.h"
#include "FastaGroup.h"
#include "Fasta.h"
#include "Fetch.h"
#include "Database.h"
#include <QVBoxLayout>
#include <QPushButton>
#include <QLabel>
#include <QLineEdit>
#include <QComboBox>
#include <iostream>

Fetch::Fetch(QWidget *parent) : QMainWindow(parent)
{
	_db = NULL;
}

void Fetch::populate()
{
	QWidget *window = new QWidget();
	QVBoxLayout *box = new QVBoxLayout();
	window->setLayout(box);

	{
		QHBoxLayout *hbox = new QHBoxLayout();
		QLabel *l = new QLabel("Country: ", window);
		hbox->addWidget(l);

		QComboBox *cb = new QComboBox();
		
		std::vector<std::string> list = _db->countryList();
		cb->addItem("");

		for (size_t i = 0; i < list.size(); i++)
		{
			QString qs = QString::fromStdString(list[i]);
			
			if (qs.length() == 0)
			{
				cb->addItem("<not specified>");
			}
			else
			{
				cb->addItem(qs);
			}
		}

		cb->setObjectName("Country");

		connect(cb, (SIGNAL(currentIndexChanged(int))), 
		        this, SLOT(updateCount()));

		hbox->addWidget(cb);

		box->addLayout(hbox);
	}
	
	{
		QHBoxLayout *hbox = new QHBoxLayout();
		QRegExp rx("[0-9]...-[0-9].-[0-9].");
		QValidator *validator = new QRegExpValidator(rx, this);

		{
			QLabel *l = new QLabel("From date: ", window);
			hbox->addWidget(l);
		}

		{
			QLineEdit *e = new QLineEdit(this);
			e->setPlaceholderText("YYYY-MM-DD");
			e->setValidator(validator);
			e->setObjectName("From");
			connect(e, &QLineEdit::editingFinished, 
			        this, &Fetch::updateCount);
			hbox->addWidget(e);
		}

		{
			QLabel *l = new QLabel("to: ", window);
			hbox->addWidget(l);
		}

		{
			QLineEdit *e = new QLineEdit(this);
			e->setPlaceholderText("YYYY-MM-DD");
			e->setValidator(validator);
			e->setObjectName("To");
			connect(e, &QLineEdit::editingFinished, 
			        this, &Fetch::updateCount);
			hbox->addWidget(e);
		}

		{
			QLabel *l = new QLabel("(inclusive)", window);
			hbox->addWidget(l);
		}

		box->addLayout(hbox);
	}
	
	{
		QHBoxLayout *hbox = new QHBoxLayout();

		{
			QLabel *l = new QLabel("Group name: ", window);
			hbox->addWidget(l);
		}

		{
			QLineEdit *e = new QLineEdit(this);
			e->setObjectName("Name");
			hbox->addWidget(e);
		}

		box->addLayout(hbox);
	}
	
	{
		QHBoxLayout *hbox = new QHBoxLayout();

		{
			QLabel *l = new QLabel("", window);
			l->setObjectName("Report");
			hbox->addWidget(l);
		}

		box->addLayout(hbox);
		
	}

	{
		QHBoxLayout *hbox = new QHBoxLayout();
		QPushButton *b = new QPushButton("Load", window);
		connect(b, &QPushButton::clicked, this, &Fetch::getSequences);
		hbox->addWidget(b);
		box->addLayout(hbox);
	}

	setCentralWidget(window);
	updateCount();
}

std::string Fetch::makeQuery(bool how_many)
{
	std::string country, from, to;

	QString qCountry = findChild<QComboBox *>("Country")->currentText();
	country = qCountry.toStdString();
	
	{
		QLineEdit *e = findChild<QLineEdit *>("From");
		QString qFrom = e->text();

		if (e->hasAcceptableInput())
		{
			from = qFrom.toStdString();
		}
		else if (qFrom.length() > 0)
		{
			std::cout << "Problem with from date" << std::endl;
			return "";
		}
	}

	{
		QLineEdit *e = findChild<QLineEdit *>("To");
		QString qTo = e->text();

		if (e->hasAcceptableInput())
		{
			to = qTo.toStdString();
		}
		else if (qTo.length() > 0)
		{
			std::cout << "Problem with to date" << std::endl;
			return "";
		}
	}
	
	bool where = false;
	
	where = (country.length() > 0) || (to.length() > 0) || (from.length() > 0);

	std::string q = "SELECT ";
	
	if (how_many)
	{
		q += "COUNT(*) AS count ";
	}
	else
	{
		q += "*, (JULIANDAY(sample_date) - JULIANDAY('2020-01-01')) ";
		q += "AS epi_days ";

	}

	q += "FROM sequences ";
	
	if (where)
	{
		q += "WHERE ";
	}
	
	if (country.length())
	{
		q += "country = '" + country + "' AND ";
	}
	
	if (from.length())
	{
		q += "sample_date >= DATE('" + from + "') AND ";
	}
	
	if (to.length())
	{
		q += "sample_date <= DATE('" + to + "') AND ";
	}
	
	for (size_t i = 0; i < 4 && where; i++)
	{
		q.pop_back();
	}

	q += ";";

	return q;
}

void Fetch::updateCount()
{
	std::string q = makeQuery(true);
	std::cout << q << std::endl;

	if (_db == NULL)
	{
		return;
	}

	_db->openConnection();
	_db->query(q);
	
	std::vector<SeqResult> results = _db->results();
	
	if (results.size())
	{
		QLabel *l = findChild<QLabel *>("Report");
		QString text = QString::fromStdString(results[0]["count"]);
		text += " sequences in database";
		l->setText(text);
	}

	_db->closeConnection();
}

void Fetch::getSequences()
{
	std::string q = makeQuery(false);

	_db->openConnection();

	_db->query(q);
	
	std::vector<SeqResult> results = _db->results();
	processResults(results);

	_db->closeConnection();
}

void Fetch::processResults(std::vector<SeqResult> results)
{
	std::string grpname = findChild<QLineEdit *>("Name")->text().toStdString();
	
	if (grpname.length() == 0)
	{
		grpname = "Untitled";
	}

	FastaGroup *group = new FastaGroup(FastaMaster::master());
	group->setCustomName(grpname);

	if (FastaMaster::master()->fastaCount() > 0)
	{
		group->addFasta(FastaMaster::master()->fasta(0));
	}

	for (size_t i = 0; i < results.size(); i++)
	{
		Fasta *f = Fasta::fastaFromDatabase(results[i]);
		group->addFasta(f);
		FastaMaster::master()->addFasta(f);
	}
	
	group->updateText();
	FastaMaster::master()->addTopLevelItem(group);
	FastaMaster::master()->setCurrentItem(group);
	group->highlightRange();
}
