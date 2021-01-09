#ifndef SPLICESITEMAP_H
#define SPLICESITEMAP_H

#include <map>
#include "bntseq.h"
//data types
enum {SPLICE_DONOR, SPLICE_ACCEPTOR, SPLICE_NONCANONICAL};

typedef struct __splice_site_t {
	uint32_t type:1, strand:1; //0: 5' or 1:3' splice site
	uint64_t pos;
	int n_partners, max_partners;
	struct __splice_site_t **partners;
} splice_site_t;

typedef std::map<uint64_t, splice_site_t *> splice_site_map_t;
typedef splice_site_map_t::iterator splice_site_map_iterator_t;


//declaration of subroutines
#ifdef __cplusplus
extern "C" {
#endif


splice_site_map_t* build_splice_site_map (const bntseq_t *bns, const char *file);
splice_site_t* retrieve_splice_site (splice_site_map_t* splice_site_map, uint64_t pos, uint32_t strand, uint32_t type);
//splice_site_t* partner_splice_site (const splice_site_t *s);
void destroy_splice_site_map (splice_site_map_t *splice_site_map_t);


#ifdef __cplusplus
}
#endif

#endif
