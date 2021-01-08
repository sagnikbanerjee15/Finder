#include <stdlib.h>
#include <string.h>

#include <iostream>
#include <map>
#include <string>

#include "BedFile.h"
#include "utils.h"
#include "splicesitemap.h"
//using namespace __gnu_cxx;

using namespace std;

//encapsulated c++ data types


//encapsulated c++ global variable
//static splice_site_map_t splice_site_map;


static uint64_t generate_splice_site_key (uint64_t pos, uint32_t strand, uint32_t type)
{
	return (pos<<2) + (strand << 1) + type;
}

static bool exist_splice_site (splice_site_map_t* splice_site_map, uint64_t k /*pos, uint32_t strand, uint32_t type*/)
{
	//uint64_t k = generate_splice_site_key (pos, strand, type);
	splice_site_map_iterator_t iter = splice_site_map->find (k);
	if (iter != splice_site_map->end()) return true;
	else return false;
}

splice_site_t* retrieve_splice_site (splice_site_map_t* splice_site_map, uint64_t pos, uint32_t strand, uint32_t type)
{
	if (splice_site_map == 0) return 0;
	uint64_t k = generate_splice_site_key (pos, strand, type);
	//splice_site_map_iterator_t iter = splice_site_map->find (k);
	if (exist_splice_site (splice_site_map, k)) return (*splice_site_map)[k];
	else return 0;
}


splice_site_t** get_splice_site_partners (const splice_site_t *s)
{
	return s->partners;
}

static void add_splice_site_partner (splice_site_t *s, splice_site_t *p)
{
	xassert (s->type != p->type, "the type of a splice site cannot be the same as its partner\n");
	if (s->n_partners >= s->max_partners) {
		s->max_partners = s->max_partners == 0 ? 1 : s->max_partners << 1;
		s->partners = (splice_site_t **) realloc (s->partners, s->max_partners * sizeof (splice_site_t *));
	}
	s->partners[s->n_partners++] = p;
}

static void add_intron (splice_site_map_t *splice_site_map, const bntann1_t *p, BedLine *intron)
{
	//int i = 0;
/*
	if ((strcmp (intron->GetChrom ().c_str(), "chr7") == 0 && intron->GetChromStart() == 135589784 && intron->GetChromEnd ()== 135591633)
		||(strcmp (intron->GetChrom ().c_str(), "chr7") == 0 && intron->GetChromStart() == 135589784 && intron->GetChromEnd ()== 135590243))
	{
		fprintf(stderr, "OK\n");
	}
*/
	uint64_t offset = p->offset;
	uint32_t strand = intron->GetStrand () == '+' ? 0 : 1; //1 for "-", and 0 for "+"
	xassert (intron->GetChromStart () < p->len, "coordinates of junctions exceeds the length of chromosome.\n");

	splice_site_t *donor, *acceptor;
	//upstream splice site
	//uint64_t pos = offset + intron->GetChromStart ();
	uint64_t k = generate_splice_site_key (offset + intron->GetChromStart (),
			strand, strand ? SPLICE_ACCEPTOR : SPLICE_DONOR);

	if (exist_splice_site (splice_site_map, k)) donor = (*splice_site_map)[k];
	else {
		donor = (splice_site_t *) calloc(1, sizeof (splice_site_t));
		donor->n_partners = donor->max_partners = 0;
		donor->strand = strand;
		donor->pos = offset + intron->GetChromStart ();
		donor->type = strand ? SPLICE_ACCEPTOR : SPLICE_DONOR;
		(*splice_site_map) [k] = donor;
		//i++;
	}


	//downstream splice site
	k = generate_splice_site_key (offset + intron->GetChromEnd (),
			strand, strand ? SPLICE_DONOR : SPLICE_ACCEPTOR);

	if (exist_splice_site (splice_site_map, k)) acceptor = (*splice_site_map)[k];
	else {
		acceptor = (splice_site_t *) calloc(1, sizeof (splice_site_t));
		acceptor->n_partners = acceptor->max_partners = 0;
		acceptor->strand = strand;
		acceptor->pos = offset + intron->GetChromEnd ();
		acceptor->type = strand ? SPLICE_DONOR : SPLICE_ACCEPTOR;
		(*splice_site_map) [k] = acceptor;
		//i++;
	}

	add_splice_site_partner (donor, acceptor);
	add_splice_site_partner (acceptor, donor);

	//return i;
}


static int splice_site_comp (const void *a, const void *b)
{
	splice_site_t **sa = (splice_site_t**) a;
	splice_site_t **sb = (splice_site_t**) b;

	return (*sa)->pos > (*sb)->pos ? 1 : ((*sa)->pos < (*sb)->pos ? -1 : 0);
}


static void sort_splice_sites (splice_site_t ** ss, int n_ss)
{
	qsort (ss, n_ss, sizeof (splice_site_t *), splice_site_comp);
}


/*
 * build a hash table of splice sites from a BED file of introns
 * so that they can be easily accessed using coordinate, strand and type
 *
 * return: the new list of splice sites
 *
 * NOTE: that alternative introns are allowed (which means duplicate splice sites),
 * but exact duplicate introns should be removed from the input
 */

struct cmp_str
{
   bool operator()(char const *a, char const *b) {
      return strcmp(a, b) < 0;
   }
};

splice_site_map_t* build_splice_site_map (const bntseq_t *bns, const char *file)
{
	//int i = 0;

	splice_site_map_t *splice_site_map = new splice_site_map_t;
	//uint64_t pacseq_pos, k;
	BedFile in (file);

	//build a hash map of chromosomes using chrom as key
	map<const char*, bntann1_t *, cmp_str> bntseq_map; //hash of chromosomes
	for (int i = 0; i != bns->n_seqs; ++i) {
		bntann1_t *p = bns->anns + i;
		bntseq_map[p->name] = p;
	}

	//add splice sites
	//i = 0;
	while (BedLine* bl = in.NextBedLine()) {
		xassert (bl->GetColNum () >= 3, "at least three columns for each intron\n");
		if (bl->GetColNum () >= 10)
			xassert (bl->GetBlockCount () == 1, "each intron can have only one block\n");

		const char * chrom = bl->GetChrom().c_str ();
		map<const char*, bntann1_t *, cmp_str>::iterator iter = bntseq_map.find (chrom);
		if (iter != bntseq_map.end()) {
			bntann1_t *p = bntseq_map[chrom];
			add_intron (splice_site_map, p, bl);
		}
		delete bl;
	}

	//sort splice site partners according to coordinates
	splice_site_map_iterator_t iter;
	for (iter = splice_site_map->begin (); iter != splice_site_map->end (); ++iter) {
		splice_site_t *s = iter->second;
		sort_splice_sites (s->partners, s->n_partners);
	}

	return splice_site_map;
}

void destroy_splice_site_map (splice_site_map_t* splice_site_map)
{

	splice_site_map_iterator_t iter;
	for (iter = splice_site_map->begin (); iter != splice_site_map->end (); ++iter) {
		splice_site_t *s = iter->second;
		free (s->partners);
		free (s);
	}
	delete splice_site_map;
}




//using namespace __gnu_cxx;
/*
int main(){
	//hash_map<const char*, int, hash<const char*>, eqstr> months;
 
	splice_site_map_t splice_site_map;
	splice_site_t *ss = (splice_site_t*) calloc(1, sizeof(splice_site_t));

	ss->type = 1;
	ss->pos = 10;
	ss->partner = 0;

	splice_site_map[ss->pos] = ss;
	
	splice_site_map_iterator_t iter = splice_site_map.find (10);

  	if (iter != splice_site_map.end()) {
  		cout << "type: " << iter->second->type << endl;
		cout << "pos:  " << iter->second->pos << endl;
	}
	free (ss);
  //cout << "april     -> " << months[4] << endl;
  //cout << "june      -> " << months[6] << endl;
  //cout << "november  -> " << months[11] << endl; 
  
  return 0;
}
*/
