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

#include "Fasta.h"
#include "FastaGroup.h"
#include "FastaMaster.h"

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <QMenu>
#include <hcsrc/Blast.h>
#include <hcsrc/FileReader.h>

bool Fasta::_justify = true;

inline bool isAddition(std::string m)
{
	return (m.back() == '+');
}

inline bool isDeletion(std::string m)
{
	return (m.back() == '-');
}

inline int residue(std::string m)
{
	int mut = atoi(&m.c_str()[1]);
	return mut;
}

inline bool isResidue(std::string m, int i)
{
	int mut = atoi(&m.c_str()[1]);
	return mut == i;
}

void Fasta::decrementResidue(std::string &m, int go_back)
{
//	unsigned char orig = m[0];
	int mut = atoi(&m.c_str()[1]);
	unsigned char last = m.back();
	std::string new_mut;
	unsigned char replacement = _ref[mut - go_back];
	new_mut += replacement;
	new_mut += i_to_str(mut - go_back);
	new_mut += last;
	m = new_mut;
}

Fasta::Fasta(std::string name)
{
	_isRef = false;
	_problematic = false;
	_compared = false;
	_name = name;
	_orf = -1;
	_stop = -1;
	_offset = -1;
	setText(0, QString::fromStdString(name));
}

unsigned char Fasta::letter(int i)
{
	return _result[_refToMe[i]];
}

void Fasta::setSequence(std::string seq, bool protein)
{
	if (protein)
	{
		_result = seq;
		return;
	}

	_seq = seq;
}

void Fasta::description()
{
	std::cout << "Sequence of " << _seq.size() << " nucleotides." << std::endl;
}

bool Fasta::findNextORF()
{
	bool success = true;
	if (true || _orf < 0)
	{
		success = findNextATG();
	}
	
	if (!success)
	{
		return false;
	}

	success = findNextStopCodon();
	
	if (!success)
	{
//		return findNextATG();
	}
	
	return success;
}

bool Fasta::findNextATG()
{
	_orf++;
	bool success = false;

	while (_orf < (int)_seq.length() - 3)
	{
		if (_seq.substr(_orf, 3) == "ATG")
		{
			_stop = _seq.length() - 1;
			success = true;
			break;
		}
		
		_orf++;
	}
	
	if (success)
	{
//		std::cout << "ATG at index: "  << _orf << std::endl;
	}
	
	return (success);
}

bool Fasta::findNextStopCodon()
{
	_stop--;
	int remainder = ((_stop - _orf) % 3);
	_stop -= remainder;
	bool success = false;

	while (_stop > _orf + 3)
	{
		std::string tmp = _seq.substr(_stop, 3);
		if (tmp == "TAA" || tmp == "TAG" || tmp == "TGA")
		{
			success = true;
			break;
		}
		
		_stop -= 3;
	}
	
	return success;
}

std::string aaLetter(std::string str)
{
	if (str == "TTT" || str == "TTC")
	{
		return "F";
	}
	else if (str == "TTA" || str == "TTG" || str == "CTT" ||
	         str == "CTC" || str == "CTA" || str == "CTG")
	{
		return "L";
	}
	else if (str == "ATT" || str == "ATC" || str == "ATA")
	{
		return "I";
	}
	else if (str == "ATG")
	{
		return "M";
	}
	else if (str == "GTT" || str == "GTC" || str == "GTA"
	         || str == "GTG")
	{
		return "V";
	}
	else if (str == "TCT" || str == "TCC" || str == "TCA"
	         || str == "TCG" || str == "AGT" || str == "AGC")
	{
		return "S";
	}
	else if (str == "CCT" || str == "CCC" || str == "CCA"
	         || str == "CCG")
	{
		return "P";
	}
	else if (str == "ACT" || str == "ACC" || str == "ACA"
	         || str == "ACG")
	{
		return "T";
	}
	else if (str == "GCT" || str == "GCC" || str == "GCA"
	         || str == "GCG")
	{
		return "A";
	}
	else if (str == "TAT" || str == "TAC")
	{
		return "Y";
	}
	else if (str == "CAT" || str == "CAC")
	{
		return "H";
	}
	else if (str == "CAA" || str == "CAG")
	{
		return "Q";
	}
	else if (str == "AAT" || str == "AAC")
	{
		return "N";
	}
	else if (str == "AAA" || str == "AAG")
	{
		return "K";
	}
	else if (str == "GAT" || str == "GAC")
	{
		return "D";
	}
	else if (str == "GAA" || str == "GAG")
	{
		return "E";
	}
	else if (str == "TGT" || str == "TGC")
	{
		return "C";
	}
	else if (str == "TGG")
	{
		return "W";
	}
	else if (str == "AGA" || str == "AGG" || str == "CGT"
	         || str == "CGC" || str == "CGA" || str == "CGG")
	{
		return "R";
	}
	else if (str == "GGA" || str == "GGG" || str == "GGT"
	         || str == "GGC")
	{
		return "G";
	}

	return " ";
}

std::string Fasta::generateSequence()
{
	std::string aaSeq = "";
	for (int i = _orf; i < _stop + 3; i += 3)
	{
		std::string tmp = _seq.substr(i, 3);
		
		std::string aa = aaLetter(tmp);
		aaSeq += aa;
	}
	
	return aaSeq;
}

bool Fasta::roughlyAlign(std::string mine, std::string ref, int minRes)
{
	if (mine == "" && hasResult())
	{
		mine = result();
	}
	else if (mine == "")
	{
		return false;
	}
	
	int size = 10;

	for (size_t i = 0; i < ref.length(); i += size / 2)
	{
		std::string probe = ref.substr(i, size);

		if ((int)probe.length() < size)
		{
			break;
		}

		if (probe.find(" ") != std::string::npos)
		{
			continue;
		}

		/* position of probe in my sequence */
		size_t loc = mine.find(probe);

		if (loc != std::string::npos)
		{
			/* how much further forward we are in relation to ref */
			int diff = (int)loc - (int)i;
			
			/* which is also the beginning...*/
			int beginning = diff;
			
			/* but we have minRes residues before then, which we may want */
			beginning -= minRes;
			
			/* but this might be too much */
			if (beginning < 0)
			{
				minRes -= beginning;
				beginning = 0;
			}
			
			/* we do want the full length of the reference */
			int end = ref.length();
			
			/* but this might be too much */
			
			if (end - beginning > (int)mine.length())
			{
				end = mine.length() - beginning;
			}

			_result = mine.substr(beginning, end);
			_offset = minRes;
			_seq.clear();
			return true;
		}
	}
	
	return false;
}

std::string Fasta::roughCompare(std::string ref, int minRes)
{
	if (hasResult())
	{
		if (roughlyAlign(_result, ref, minRes))
		{
			return _result;
		}

		return _result;
	}

	int count = 0;
	while (findNextORF())
	{
		count++;
		std::string seq = generateSequence();
		
		if (roughlyAlign(seq, ref, minRes))
		{
			std::cout << "After trying " << count << " ORFs, "
			" length " << seq.length() <<  " nt..." << std::endl;
			return seq;
		}
	}

	std::cout << "Did not find" << std::endl;
	return "";
}

int Fasta::gapCount()
{
	int count = 0;
	char ch = ' ';
	for (size_t i = 0; (i = _result.find(ch, i)) != std::string::npos; i++) 
	{
		count++;
	}

	return count;
}

void Fasta::addMutation(char fromWhat, int mut, char towhat)
{
	std::ostringstream ss;
	ss << fromWhat << mut + _offset << towhat;
	std::string str = ss.str();

	_mutations.push_back(str);
}

void Fasta::clearMutations()
{
	_mutations.clear();
	_compared = false;
}

void Fasta::carefulCompareWithFasta(Fasta *f)
{
	if (!hasResult() && !f->hasResult())
	{
		return;
	}
	
	if (f == this)
	{
		_isRef = true;
		_ref = _result;
		_offset = 0;
		organiseMap();
		findGlycosylations();
		return;
	}

	carefulCompareWithString(f->result());
	organiseMap();
	findGlycosylations();
	removeDuplicateGlycosylations(f);

	for (size_t j = 0; j < mutationCount(); j++)
	{
		if (mutation(j).back() == '>' || mutation(j).back() == '<')
		{
			std::cout << name() << " " << mutation(j) << std::endl;
		}
	}
}

void Fasta::removeDuplicateGlycosylations(Fasta *f)
{
	for (size_t j = 0; j < f->mutationCount(); j++)
	{
		std::string refMut = f->mutation(j);
		bool found = false;
		for (size_t i = 0; i < mutationCount(); i++)
		{
			std::string myMut = mutation(i);
			
			if (myMut.back() == '<')
			{
				continue;
			}

			if (myMut == refMut)
			{
				found = true;
				_mutations.erase(_mutations.begin() + i);
				break;
			}
		}
		
		if (!found && refMut.back() == '>')
		{
			refMut.back() = '<';
			_mutations.push_back(refMut);
		}
	}

}

void Fasta::writeAlignment(std::ofstream &file)
{
	file << ">" << name() << std::endl;
	file << _left << std::endl;
	file << _align << std::endl;
	file << _right << std::endl;
	file << std::endl;
}

bool trio_can_glycosylate(std::string trio)
{
	bool ok = true;
	
	if (trio[0] != 'N' || trio[1] == 'P')
	{
		ok = false;
	}
	if (trio[2] != 'S' && trio[2] != 'T')
	{
		ok = false;
	}
	
	return ok;
}

void Fasta::findGlycosylations()
{
	return;
	for (size_t mut = 0; mut < _ref.length() - 3; mut++)
	{
		size_t i = _refToMe[mut];
		std::string trio = "   ";
		
		if (i < _result.length())
		{
			trio = _result.substr(i, 3);
		}
		
		/* if we have white spaces, we have to make something compatible */
		if (!_isRef && trio.find(' ') != std::string::npos)
		{
			std::string orig = _ref.substr(mut, 3);

			if (trio_can_glycosylate(orig))
			{
				if (trio[0] == ' ') trio[0] = 'N';
				if (trio[2] == ' ') trio[2] = 'T';
			}
		}

		if (trio_can_glycosylate(trio))
		{
			addMutation('N', mut, '>');
		}
	}
}

void Fasta::carefulCompareWithString(std::string seq2)
{
	std::string seq1 = result();
	
	_mutations.clear();
	
	if (seq1.length() == 0 || seq2.length() == 0)
	{
		return;
	}

	_ref = seq2;

	int muts, dels;
	Alignment ala, alb;
	setup_alignment(&ala, "");
	setup_alignment(&alb, "");
	compare_sequences_and_alignments(seq1, seq2, &muts, &dels, ala, alb, 2);
	tidy_alignments(ala, alb);

	std::ostringstream ssleft, ssalign, ssright;
	std::vector<int> indices;
	print_alignments(ala, alb, ssleft, ssalign, ssright, indices);
	
	delete_alignment(&ala);
	delete_alignment(&alb);
	
	_left =   ssleft.str();
	_align = ssalign.str();
	_right = ssright.str();
	
	/*
	std::cout << _left << std::endl;
	std::cout << _align << std::endl;
	std::cout << _right << std::endl;
	*/
	
	int running = 0;
	int offset = 0;
	int plus = 0;
	bool ending = false;

	for (size_t i = 0; i < _align.size(); i++)
	{
		if (_align[i] == '.' || _align[i] == '*')
		{
			offset += running;
			running = 0;
		}

		offset = 0;

		if (_align[i] == '.')
		{
			continue;
		}

		if (_align[i] == '*')
		{
			addMutation(_right[i], indices[i] + offset, _left[i]);
			continue;
		}
		
		if (_align[i] == '+' && !ending)
		{
			addMutation(_left[i], indices[i] + offset, '+');
			running--;
			plus++;
		}
		
		if (_align[i] == '-')
		{
			addMutation(_right[i], indices[i] + offset, '-');
			running++;
		}
		
		if (i >= _align.size() - 10)
		{
			ending = true;
		}
		
		if (plus > 10)
		{
			_problematic = true;
		}
	}
	
	leftJustifyDeletions();
	
	_compared = true;
	
	FastaMaster::master()->addValue(this, "mutations", mutationSummary());
}

void Fasta::leftJustifyDeletions()
{
	if (!_justify)
	{
		return;
	}

	sortMutations();
	
	for (size_t i = 0; i < mutationCount(); i++)
	{
		std::string m = mutation(i);
		int start = residue(m);
		
		if (!isDeletion(mutation(i)))
		{
			continue;
		}
		
		int end = start + 1;
		int count = 1;

		for (size_t j = i + 1; j < mutationCount(); j++)
		{
			if (!isResidue(mutation(j), end))
			{
				break;
			}
			count++; end++;
		}
		
		if (start == _offset)
		{
			continue;
		}
		
		if (end == _offset + (int)result().length() - 1)
		{
			continue;
		}
		
		std::string comparison = deletionSequence(start, end, 0);
		
		if (comparison.length() == 0)
		{
			continue;
		}
		
		if (comparison.find(' ') != std::string::npos)
		{
			continue;
		}

		int go_back = 0;

		while (true)
		{
			go_back++;
			std::string back_one = deletionSequence(start, end, go_back);

			if (back_one.find(' ') != std::string::npos)
			{
				break;
			}

			if (comparison.length() != back_one.length() ||
			    comparison != back_one)
			{
				break;
			}
		}

		go_back--;
		
		if (go_back == 0)
		{
			i += count - 1;
			continue;
		}
		
		for (size_t j = i; j < i + count; j++)
		{
			std::string mut = _mutations[j];
			decrementResidue(mut, go_back);
			_mutations[j] = mut;
		}

		i += count - 1;
	}
}

std::string Fasta::deletionSequence(int start, int end, int go_back)
{
	int left = start - 2;
	int right = end + 2;
	std::string seq;
	
	start -= go_back;
	end -= go_back;

	for (int i = left; i <= right; i++)
	{
		if (i >= start && i < end)
		{
			continue;
		}

		if (i < 0 || i >= (int)_ref.length())
		{
			seq += "?";
			continue;
		}

		seq += _ref[i];
	}
	
	return seq;
}

std::string Fasta::mutationSummary()
{
	if (mutationCount() == 0)
	{
		return " ";
	}

	std::ostringstream ss;
	
	for (size_t i = 0; i < mutationCount(); i++)
	{
		ss << mutation(i) << " ";
	}

	std::string str = ss.str();
	
	if (str.length() > 0)
	{
		str.pop_back();
	}

	return str;
}

void Fasta::nudgeMap(int start, int dir)
{
	for (IntMap::iterator it = _refToMe.begin(); it != _refToMe.end(); it++)
	{
		if (it->first >= start)
		{
			it->second += dir;
		}
	}
}

void Fasta::organiseMap()
{
	_refToMe.clear();
	_meToRef.clear();

	for (size_t i = 0; i < _result.length(); i++)
	{
		_refToMe[i + _offset] = i;
	}

	for (size_t i = 0; i < mutationCount(); i++)
	{
		std::string m = mutation(i);
		int resi = residue(m);
		
		if (isAddition(m))
		{
			nudgeMap(resi, +1);
		}
		else if (isDeletion(m))
		{
			nudgeMap(resi, -1);
		}
	}

	for (IntMap::iterator it = _refToMe.begin(); it != _refToMe.end(); it++)
	{
		if (_meToRef.count(it->second) == 0)
		{
			_meToRef[it->second] = it->first;
		}
	}
}

void Fasta::sortMutations()
{
	std::vector<MutInt> tmp;

	for (size_t i = 0; i < mutationCount(); i++)
	{
		MutInt mi;
		mi.mut = mutation(i);
		mi.resi = residue(mi.mut);
		tmp.push_back(mi);
	}
	
	std::sort(tmp.begin(), tmp.end(), resi_less);
	_mutations.clear();

	for (size_t i = 0; i < tmp.size(); i++)
	{
		_mutations.push_back(tmp[i].mut);
	}
}

void Fasta::loadMutations(std::string muts, std::string ref)
{
	_compared = true;
	_ref = ref;

	_mutations.clear();
	
	std::vector<std::string> components = split(muts, ' ');
	int plus = 0;
	
	for (size_t i = 0; i < components.size(); i++)
	{
		if (components[i].length() == 0)
		{
			continue;
		}
		
		if (components[i].back() == '+')
		{
			plus++;
		}

		_mutations.push_back(components[i]);
	}

	if (plus > 10)
	{
		_problematic = true;
	}
	
	leftJustifyDeletions();
}

void Fasta::setIsProblematic()
{
	_problematic = true;
	emit refreshMutations();
}

void Fasta::giveMenu(QMenu *m, FastaGroup *g)
{
	QAction *act = m->addAction("Set as problematic");
	connect(act, &QAction::triggered, this, &Fasta::setIsProblematic);
}

int Fasta::oneSidedMutations(Fasta *other)
{
	int total = mutationCount();

	for (size_t i = 0; i < mutationCount(); i++)
	{
		bool found = other->hasMutation(mutation(i));
		
		if (found)
		{
			total--;
		}
	}

	return total;
}

int Fasta::sharedMutations(Fasta *other)
{
	int total = 0;
	total += oneSidedMutations(other);
	total += other->oneSidedMutations(this);

	return total;
}

bool Fasta::hasMutation(std::string mut)
{
	bool found = (std::find(_mutations.begin(), _mutations.end(), mut)
	              != _mutations.end());
	
	return found;
}

std::string Fasta::selectQuery()
{
	std::string me;
	me = "SELECT * FROM sequences WHERE name = '" + name() + "';";
	return me;
}

std::string Fasta::updateQuery()
{
	std::string me;

	me = "UPDATE sequences ";
	me += "SET ";
	me += "source = '" + _source + "', ";
	me += "country = '" + _country + "', ";
	me += "sample_date = '";
	me += FastaMaster::master()->valueForKey(this, "sample_date") + "', "; 
	me += "added_date = DATE('now'), ";
	me += "protein_sequence = '" + result() + "' ";
	
	std::string muts = mutationSummary();
	
	if (hasCompared())
	{
		me += ", mutations = '" + mutationSummary() + "' ";
	}

	me += "WHERE name = '" + name() + "'; ";

	return me;
}

std::string Fasta::insertQuery()
{
	std::string me;
	me = "INSERT OR IGNORE INTO sequences (name) VALUES ('" + name() + "');";

	/*
	me = "INSERT INTO sequences (";
	me += "name, ";
	me += "source, ";
	me += "country, ";
	me += "sample_date, ";
	me += "added_date, ";
	me += "protein_sequence, ";
	me += "mutations";
	me += ") VALUES (";
	me += "'" + name() + "', ";
	me += "'" + source() + "', ";
	me += "'" + country() + "', ";
	me += "'" + FastaMaster::master()->valueForKey(this, "sample_date") + "', ";
	me += "DATE('now'), ";
	me += "'" + result() + "', ";
	me += "'" + mutationSummary() + "');";
	*/

	return me;
}

void Fasta::setCountry(std::string country)
{
	for (size_t i = 0; i < country.size(); i++)
	{
		if (country[i] == ' ')
		{
			country.erase(country.begin() + i);
			i--;
		}
	}
	
	to_lower(country);
	
	_country = country;
}

void Fasta::figureOutFromName()
{
	std::string n = name();
	std::vector<std::string> bits = split(name(), '/');
	
	if (_source == "cog" && bits.size() >= 1)
	{
		_country = bits[0];
	}
	else if (_source == "gisaid" && bits.size() >= 2)
	{
		_country = bits[1];
	}
	
	bits = split(name(), '|');
	
	if (_source == "gisaid" && bits.size() >= 3)
	{
		std::string value = bits[2];
		FastaMaster::master()->addValue(this, "sample_date", value);
	}
}

Fasta *Fasta::fastaFromDatabase(SeqResult &r)
{
	Fasta *f = new Fasta(r["name"]);
	f->setSource(r["source"]);
	f->setCountry(r["country"]);
	f->setSequence(r["protein_sequence"], true);

	std::ostringstream str;
	str << std::setfill('0') << std::setw(7) << r["epi_days"];

	FastaMaster::master()->addValue(f, "sample_date", r["sample_date"]);
	FastaMaster::master()->addValue(f, "mutations", r["mutations"]);
	FastaMaster::master()->addValue(f, "epi_days", str.str());
	
	if (r["mutations"].length() > 1)
	{
		f->loadMutations(r["mutations"], 
		                 FastaMaster::master()->fasta(0)->result());
	}

	return f;
}

void Fasta::fetchValueForTitle(std::string title)
{
	_lastValue = FastaMaster::master()->valueForKey(this, title);
}
