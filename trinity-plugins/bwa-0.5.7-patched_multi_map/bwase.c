#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "stdaln.h"
#include "bwase.h"
#include "bwtaln.h"
#include "bntseq.h"
#include "utils.h"
#include "kstring.h"

static int g_log_n[256];

void bwa_aln2seq_core(int n_aln, const bwt_aln1_t *aln, bwa_seq_t *s, int set_main, int n_multi)
{
	int i, cnt, best;
	if (n_aln == 0) {
		s->type = BWA_TYPE_NO_MATCH;
		s->c1 = s->c2 = 0;
		return;
	}

	if (set_main) {
		best = aln[0].score;
		for (i = cnt = 0; i < n_aln; ++i) {
			const bwt_aln1_t *p = aln + i;
			if (p->score > best) break;
			if (drand48() * (p->l - p->k + 1) > (double)cnt) {
				s->n_mm = p->n_mm; s->n_gapo = p->n_gapo; s->n_gape = p->n_gape; s->strand = p->a;
				s->score = p->score;
				s->sa = p->k + (bwtint_t)((p->l - p->k + 1) * drand48());
			}
			cnt += p->l - p->k + 1;
		}
		s->c1 = cnt;
		for (; i < n_aln; ++i) cnt += aln[i].l - aln[i].k + 1;
		s->c2 = cnt - s->c1;
		s->type = s->c1 > 1? BWA_TYPE_REPEAT : BWA_TYPE_UNIQUE;
	}

	if (n_multi) {
		int k, rest, n_occ, z = 0;
		for (k = n_occ = 0; k < n_aln; ++k) {
			const bwt_aln1_t *q = aln + k;
			n_occ += q->l - q->k + 1;
		}
		if (s->multi) free(s->multi);
		if (n_occ > n_multi + 1) { // if there are too many hits, generate none of them
			s->multi = 0; s->n_multi = 0;
			return;
		}
		/* The following code is more flexible than what is required
		 * here. In principle, due to the requirement above, we can
		 * simply output all hits, but the following samples "rest"
		 * number of random hits. */
		rest = n_occ > n_multi + 1? n_multi + 1 : n_occ; // find one additional for ->sa
		s->multi = calloc(rest, rest * sizeof(bwt_multi1_t));
		for (k = 0; k < n_aln; ++k) {
			const bwt_aln1_t *q = aln + k;
			if (q->l - q->k + 1 <= rest) {
				bwtint_t l;
				for (l = q->k; l <= q->l; ++l) {
					s->multi[z].pos = l;
					s->multi[z].gap = q->n_gapo + q->n_gape;
					s->multi[z].mm = q->n_mm;
					s->multi[z++].strand = q->a;
				}
				rest -= q->l - q->k + 1;
			} else { // Random sampling (http://code.activestate.com/recipes/272884/). In fact, we never come here. 
				int j, i, k;
				for (j = rest, i = q->l - q->k + 1, k = 0; j > 0; --j) {
					double p = 1.0, x = drand48();
					while (x < p) p -= p * j / (i--);
					s->multi[z].pos = q->l - i;
					s->multi[z].gap = q->n_gapo + q->n_gape;
					s->multi[z].mm = q->n_mm;
					s->multi[z++].strand = q->a;
				}
				rest = 0;
				break;
			}
		}
		s->n_multi = z;
		for (k = z = 0; k < s->n_multi; ++k)
			if (s->multi[k].pos != s->sa)
				s->multi[z++] = s->multi[k];
		s->n_multi = z < n_multi? z : n_multi;
	}
}

void bwa_aln2seq(int n_aln, const bwt_aln1_t *aln, bwa_seq_t *s)
{
	bwa_aln2seq_core(n_aln, aln, s, 1, 0);
}

/* *aln points to alignments found for the current sequence, n_aln is the 
   size of the array pointed to by *aln. Array *s of size n_seq is a collection
   of SAM records that must be replicas, i.e. initialized with the same current sequence.
   This method updates sequence records in *s with placements recorded in *aln,
   until all n_seq replicas are updated. For each separate alignment record in array *aln,
   each placement corresponding to this record will be assigned to a separate record in *s 
   if there are enough elements remaining in *s, otherwise a random subset of the placements will
   be assigned to the remaining elements in *s. The total number of best placements and total
   number of non-best placements will be computed from the whole array *aln (regardless of whether it
   fits completely into *s or not) and assigned to each updated record in *s.
*/
void bwa_aln2seq_all(int n_aln, const bwt_aln1_t *aln, int n_seq, bwa_seq_t *s)
{
  int i, cnt1, cnt2, j, best, N;
	if (n_aln == 0) {
          /* there is no match found for *s  */
		s->type = BWA_TYPE_NO_MATCH;
		s->c1 = s->c2 = 0;
		return;
	}
        
        N = n_seq; // remember the size of the array
        best = aln[0].score;

        cnt1 = 0; // total number of already processed alignments (i.e. distinct placements, NOT alignment records) with best score
        cnt2 = 0; // total number of already processed alignments with inferior score(s)
	for (i = 0; i < n_aln && n_seq > 0 ; ++i) {
          const bwt_aln1_t *p_aln = aln + i;

          int N_aligns = p_aln->l-p_aln->k +1 ; // number of placements (alignments) in the current alignment record p_aln

          if (N_aligns <= n_seq) { /* we have space to save all the alignments stored in 'p_aln' */

              for ( j = 0 ; j < N_aligns ; j++ ) {
            
                  bwa_seq_t * seq = s + cnt1+ cnt2+j ;
                  seq->n_mm = p_aln->n_mm; seq->n_gapo = p_aln->n_gapo; seq->n_gape = p_aln->n_gape; seq->strand = p_aln->a;
                  seq->score = p_aln->score;
                  seq->sa = p_aln->k + j;
              }
              n_seq -= N_aligns; // we have n_seq slots remaining to store more alignments
          } else { // See also: http://code.activestate.com/recipes/272884/
            // we have to truncate, so let's select few remaining alignments randomly:
              int xj, xi, xk;
              for (xj = n_seq, xi = N_aligns, xk = 0; xj > 0; --xj, ++xk) {
                double p = 1.0, x = drand48();
                while (x < p) p -= p * xj / (xi--);

                bwa_seq_t * seq = s+cnt1+cnt2+xk ;
                seq->n_mm = p_aln->n_mm; seq->n_gapo = p_aln->n_gapo; seq->n_gape = p_aln->n_gape; seq->strand = p_aln->a;
                seq->score = p_aln->score;
                seq->sa = p_aln->l - xi;

              }
              n_seq = 0;
          }
          // cnt1 + cnt2 is the total count of hits processed so far:
          if ( p_aln->score == best  ) cnt1 += N_aligns; // we found N_aligns more placements with best score
          else cnt2 += N_aligns; // N_aligns more placements with inferior score
        }
          
        // we filled all available slots in the array *s, but there can be more alignments
        // left; we need to count them:
	for (; i < n_aln; ++i) {
          if ( aln[i].score == best  ) cnt1 += aln[i].l-aln[i].k+1;
          else cnt2 += aln[i].l-aln[i].k+1;
        }

        // now cnt1 is the total number of found alignments (placements) with best score
        // and cnt2 is the total number of found placements with worse score

        /* set counts and flags for all hits: */
	for (i = 0; i < N ; ++i) {
          bwa_seq_t * seq = s+i ;

          seq->c1 = cnt1;
          seq->c2 = cnt2;
          seq->type = seq->c1 > 1? BWA_TYPE_REPEAT : BWA_TYPE_UNIQUE;
        }
}



int bwa_approx_mapQ(const bwa_seq_t *p, int mm)
{
	int n;
	if (p->c1 == 0) return 23;
	if (p->c1 > 1) return 0;
	if (p->n_mm == mm) return 25;
	if (p->c2 == 0) return 37;
	n = (p->c2 >= 255)? 255 : p->c2;
	return (23 < g_log_n[n])? 0 : 23 - g_log_n[n];
}

/**
 * Derive the actual position in the read from the given suffix array
 * coordinates. Note that the position will be approximate based on
 * whether indels appear in the read and whether calculations are
 * performed from the start or end of the read.
 */
void bwa_cal_pac_pos_core(const bwt_t *forward_bwt, const bwt_t *reverse_bwt, int n_seqs, bwa_seq_t *s, const int max_mm, const float fnr)
{
	int max_diff;
        bwa_seq_t *seq;
        int i;

        for ( i = 0 ; i < n_seqs ; i++ ) {
          seq = s + i;
          if (seq->type != BWA_TYPE_UNIQUE && seq->type != BWA_TYPE_REPEAT) continue;

          max_diff = fnr > 0.0? bwa_cal_maxdiff(seq->len, BWA_AVG_ERR, fnr) : max_mm;
          if (seq->strand) { // reverse strand only
            seq->pos = bwt_sa(forward_bwt, seq->sa);
            seq->seQ = seq->mapQ = bwa_approx_mapQ(seq, max_diff);
          } else { // forward strand only
            /* NB: For gapped alignment, p->pos may not be correct, which
             *     will be fixed in refine_gapped_core(). This line also
             *     determines the way "x" is calculated in
             *     refine_gapped_core() when (ext < 0 && is_end == 0). */
            seq->pos = reverse_bwt->seq_len - (bwt_sa(reverse_bwt, seq->sa) + seq->len);
            seq->seQ = seq->mapQ = bwa_approx_mapQ(seq, max_diff);
          }
        }
}

void bwa_cal_pac_pos(const char *prefix, int n_seqs, bwa_seq_t *seqs, int max_mm, float fnr)
{
	int i, j;
	char str[1024];
	bwt_t *bwt;
	// load forward SA
	strcpy(str, prefix); strcat(str, ".bwt");  bwt = bwt_restore_bwt(str);
	strcpy(str, prefix); strcat(str, ".sa"); bwt_restore_sa(str, bwt);
	for (i = 0; i != n_seqs; ++i) {
          if (seqs[i].strand) bwa_cal_pac_pos_core(bwt, 0, 1, &seqs[i], max_mm, fnr);
		for (j = 0; j < seqs[i].n_multi; ++j) {
			bwt_multi1_t *p = seqs[i].multi + j;
			if (p->strand) p->pos = bwt_sa(bwt, p->pos);
		}
	}
	bwt_destroy(bwt);
	// load reverse BWT and SA
	strcpy(str, prefix); strcat(str, ".rbwt"); bwt = bwt_restore_bwt(str);
	strcpy(str, prefix); strcat(str, ".rsa"); bwt_restore_sa(str, bwt);
	for (i = 0; i != n_seqs; ++i) {
          if (!seqs[i].strand) bwa_cal_pac_pos_core(0, bwt, 1, &seqs[i], max_mm, fnr);
		for (j = 0; j < seqs[i].n_multi; ++j) {
			bwt_multi1_t *p = seqs[i].multi + j;
			if (!p->strand) p->pos = bwt->seq_len - (bwt_sa(bwt, p->pos) + seqs[i].len);
		}
	}
	bwt_destroy(bwt);
}

/* is_end_correct == 1 if (*pos+len) gives the correct coordinate on
 * forward strand. This happens when p->pos is calculated by
 * bwa_cal_pac_pos(). is_end_correct==0 if (*pos) gives the correct
 * coordinate. This happens only for color-converted alignment. */
static bwa_cigar_t *refine_gapped_core(bwtint_t l_pac, const ubyte_t *pacseq, int len, const ubyte_t *seq, bwtint_t *_pos,
									int ext, int *n_cigar, int is_end_correct)
{
	bwa_cigar_t *cigar = 0;
	ubyte_t *ref_seq;
	int l = 0, path_len, ref_len;
	AlnParam ap = aln_param_bwa;
	path_t *path;
	int64_t k, __pos = *_pos > l_pac? (int64_t)((int32_t)*_pos) : *_pos;

	ref_len = len + abs(ext);
	if (ext > 0) {
		ref_seq = (ubyte_t*)calloc(ref_len, 1);
		for (k = __pos; k < __pos + ref_len && k < l_pac; ++k)
			ref_seq[l++] = pacseq[k>>2] >> ((~k&3)<<1) & 3;
	} else {
		int64_t x = __pos + (is_end_correct? len : ref_len);
		ref_seq = (ubyte_t*)calloc(ref_len, 1);
		for (l = 0, k = x - ref_len > 0? x - ref_len : 0; k < x && k < l_pac; ++k)
			ref_seq[l++] = pacseq[k>>2] >> ((~k&3)<<1) & 3;
	}
	path = (path_t*)calloc(l+len, sizeof(path_t));

	aln_global_core(ref_seq, l, (ubyte_t*)seq, len, &ap, path, &path_len);
	cigar = bwa_aln_path2cigar(path, path_len, n_cigar);
	
	if (ext < 0 && is_end_correct) { // fix coordinate for reads mapped on the forward strand
		for (l = k = 0; k < *n_cigar; ++k) {
			if (__cigar_op(cigar[k]) == FROM_D) l -= __cigar_len(cigar[k]);
			else if (__cigar_op(cigar[k]) == FROM_I) l += __cigar_len(cigar[k]);
		}
		__pos += l;
	}

	if (__cigar_op(cigar[0]) == FROM_D) { // deletion at the 5'-end
		__pos += __cigar_len(cigar[0]);
		for (k = 0; k < *n_cigar - 1; ++k) cigar[k] = cigar[k+1];
		--(*n_cigar);
	}
	if (__cigar_op(cigar[*n_cigar-1]) == FROM_D) --(*n_cigar); // deletion at the 3'-end

	// change "I" at either end of the read to S. just in case. This should rarely happen...
	if (__cigar_op(cigar[*n_cigar-1]) == FROM_I) cigar[*n_cigar-1] = __cigar_create(3, (__cigar_len(cigar[*n_cigar-1])));
	if (__cigar_op(cigar[0]) == FROM_I) cigar[0] = __cigar_create(3, (__cigar_len(cigar[0])));

	*_pos = (bwtint_t)__pos;
	free(ref_seq); free(path);
	return cigar;
}

char *bwa_cal_md1(int n_cigar, bwa_cigar_t *cigar, int len, bwtint_t pos, ubyte_t *seq,
				  bwtint_t l_pac, ubyte_t *pacseq, kstring_t *str, int *_nm)
{
	bwtint_t x, y;
	int z, u, c, nm = 0;
	str->l = 0; // reset
	x = pos; y = 0;
	if (cigar) {
		int k, l;
		for (k = u = 0; k < n_cigar; ++k) {
			l = __cigar_len(cigar[k]);
			if (__cigar_op(cigar[k]) == FROM_M) {
				for (z = 0; z < l && x+z < l_pac; ++z) {
					c = pacseq[(x+z)>>2] >> ((~(x+z)&3)<<1) & 3;
					if (c > 3 || seq[y+z] > 3 || c != seq[y+z]) {
						ksprintf(str, "%d", u);
						kputc("ACGTN"[c], str);
						++nm;
						u = 0;
					} else ++u;
				}
				x += l; y += l;
/*		        } else if (cigar[k]>>14 == FROM_I || cigar[k]>>14 == 3) { */
                        } else if (__cigar_op(cigar[k]) == FROM_I || __cigar_op(cigar[k]) == FROM_S) {
				y += l;
				if (__cigar_op(cigar[k]) == FROM_I) nm += l;
			} else if (__cigar_op(cigar[k]) == FROM_D) {
				ksprintf(str, "%d", u);
				kputc('^', str);
				for (z = 0; z < l && x+z < l_pac; ++z)
					kputc("ACGT"[pacseq[(x+z)>>2] >> ((~(x+z)&3)<<1) & 3], str);
				u = 0;
				x += l; nm += l;
			}
		}
	} else { // no gaps
		for (z = u = 0; z < (bwtint_t)len; ++z) {
			c = pacseq[(x+z)>>2] >> ((~(x+z)&3)<<1) & 3;
			if (c > 3 || seq[y+z] > 3 || c != seq[y+z]) {
				ksprintf(str, "%d", u);
				kputc("ACGTN"[c], str);
				++nm;
				u = 0;
			} else ++u;
		}
	}
	ksprintf(str, "%d", u);
	*_nm = nm;
	return strdup(str->s);
}

void bwa_correct_trimmed(bwa_seq_t *s)
{
	if (s->len == s->full_len) return;
	if (s->strand == 0) { // forward
		if (s->cigar && __cigar_op(s->cigar[s->n_cigar-1]) == FROM_S) { // the last is S
			s->cigar[s->n_cigar-1] += s->full_len - s->len;
		} else {
			if (s->cigar == 0) {
				s->n_cigar = 2;
				s->cigar = calloc(s->n_cigar, sizeof(bwa_cigar_t));
				s->cigar[0] = __cigar_create(0, s->len);
			} else {
				++s->n_cigar;
				s->cigar = realloc(s->cigar, s->n_cigar * sizeof(bwa_cigar_t));
			}
			s->cigar[s->n_cigar-1] = __cigar_create(3, (s->full_len - s->len));
		}
	} else { // reverse
		if (s->cigar && __cigar_op(s->cigar[0]) == FROM_S) { // the first is S
			s->cigar[0] += s->full_len - s->len;
		} else {
			if (s->cigar == 0) {
				s->n_cigar = 2;
				s->cigar = calloc(s->n_cigar, sizeof(bwa_cigar_t));
				s->cigar[1] = __cigar_create(0, s->len);
			} else {
				++s->n_cigar;
				s->cigar = realloc(s->cigar, s->n_cigar * sizeof(bwa_cigar_t));
				memmove(s->cigar + 1, s->cigar, (s->n_cigar-1) * sizeof(bwa_cigar_t));
			}
			s->cigar[0] = __cigar_create(3, (s->full_len - s->len));
		}
	}
	s->len = s->full_len;
}

void bwa_refine_gapped(const bntseq_t *bns, int n_seqs, bwa_seq_t *seqs, ubyte_t *_pacseq, bntseq_t *ntbns)
{
	ubyte_t *pacseq, *ntpac = 0;
	int i, j;
	kstring_t *str;

	if (ntbns) { // in color space
		ntpac = (ubyte_t*)calloc(ntbns->l_pac/4+1, 1);
		rewind(ntbns->fp_pac);
		fread(ntpac, 1, ntbns->l_pac/4 + 1, ntbns->fp_pac);
	}

	if (!_pacseq) {
		pacseq = (ubyte_t*)calloc(bns->l_pac/4+1, 1);
		rewind(bns->fp_pac);
		fread(pacseq, 1, bns->l_pac/4+1, bns->fp_pac);
	} else pacseq = _pacseq;
	for (i = 0; i != n_seqs; ++i) {
		bwa_seq_t *s = seqs + i;
		for (j = 0; j < s->n_multi; ++j) {
			bwt_multi1_t *q = s->multi + j;
			int n_cigar;
			if (q->gap == 0) continue;
			q->cigar = refine_gapped_core(bns->l_pac, pacseq, s->len, q->strand? s->rseq : s->seq, &q->pos,
										  (q->strand? 1 : -1) * q->gap, &n_cigar, 1);
			q->n_cigar = n_cigar;
		}
		if (s->type == BWA_TYPE_NO_MATCH || s->type == BWA_TYPE_MATESW || s->n_gapo == 0) continue;
		s->cigar = refine_gapped_core(bns->l_pac, pacseq, s->len, s->strand? s->rseq : s->seq, &s->pos,
									  (s->strand? 1 : -1) * (s->n_gapo + s->n_gape), &s->n_cigar, 1);
	}

	if (ntbns) { // in color space
		for (i = 0; i < n_seqs; ++i) {
			bwa_seq_t *s = seqs + i;
			bwa_cs2nt_core(s, bns->l_pac, ntpac);
			for (j = 0; j < s->n_multi; ++j) {
				bwt_multi1_t *q = s->multi + j;
				int n_cigar;
				if (q->gap == 0) continue;
				free(q->cigar);
				q->cigar = refine_gapped_core(bns->l_pac, ntpac, s->len, q->strand? s->rseq : s->seq, &q->pos,
											  (q->strand? 1 : -1) * q->gap, &n_cigar, 0);
				q->n_cigar = n_cigar;
			}
			if (s->type != BWA_TYPE_NO_MATCH && s->cigar) { // update cigar again
				free(s->cigar);
				s->cigar = refine_gapped_core(bns->l_pac, ntpac, s->len, s->strand? s->rseq : s->seq, &s->pos,
											  (s->strand? 1 : -1) * (s->n_gapo + s->n_gape), &s->n_cigar, 0);
			}
		}
	}

	// generate MD tag
	str = (kstring_t*)calloc(1, sizeof(kstring_t));
	for (i = 0; i != n_seqs; ++i) {
		bwa_seq_t *s = seqs + i;
		if (s->type != BWA_TYPE_NO_MATCH) {
			int nm;
			s->md = bwa_cal_md1(s->n_cigar, s->cigar, s->len, s->pos, s->strand? s->rseq : s->seq,
								bns->l_pac, ntbns? ntpac : pacseq, str, &nm);
			s->nm = nm;
		}
	}
	free(str->s); free(str);

	// correct for trimmed reads
	if (!ntbns) // trimming is only enabled for Illumina reads
		for (i = 0; i < n_seqs; ++i) bwa_correct_trimmed(seqs + i);

	if (!_pacseq) free(pacseq);
	free(ntpac);
}

int64_t pos_end(const bwa_seq_t *p)
{
	if (p->cigar) {
		int j;
		int64_t x = p->pos;
		for (j = 0; j != p->n_cigar; ++j) {
			int op = __cigar_op(p->cigar[j]);
			if (op == 0 || op == 2) x += __cigar_len(p->cigar[j]);
		}
		return x;
	} else return p->pos + p->len;
}

int64_t pos_end_multi(const bwt_multi1_t *p, int len) // analogy to pos_end()
{
	if (p->cigar) {
		int j;
		int64_t x = p->pos;
		for (j = 0; j != p->n_cigar; ++j) {
			int op = __cigar_op(p->cigar[j]);
			if (op == 0 || op == 2) x += __cigar_len(p->cigar[j]);
		}
		return x;
	} else return p->pos + len;
}

static int64_t pos_5(const bwa_seq_t *p)
{
	if (p->type != BWA_TYPE_NO_MATCH)
		return p->strand? pos_end(p) : p->pos;
	return -1;
}

/* Prints <bases>\t<quals> of the sequence *p into STDOUT; 
*/
void bwa_print_seq_and_qual(bwa_seq_t *p) {
  int j;
  ubyte_t * char_ptr;

  if (p->strand == 0)
    for (j = 0; j != p->full_len; ++j) putchar("ACGTN"[(int)p->seq[j]]);
  else for (j = 0; j != p->full_len; ++j) putchar("TGCAN"[p->seq[p->full_len - 1 - j]]);
  putchar('\t');
  if (p->qual) {
    if (p->strand) {
      //      seq_reverse(p->len, p->qual, 0); // reverse quality
      for ( char_ptr = p->qual + p->len - 1 ; char_ptr >= p->qual ; char_ptr-- ) putchar( *char_ptr );
    } else {
      printf("%s", p->qual);
    }
  } else printf("*");
  
}



void bwa_print_sam1(const bntseq_t *bns, bwa_seq_t *p, const bwa_seq_t *mate, int mode, int max_top2)
{

  //  if ( strcmp(p->name,"62A73AAXX100903:1:100:10014:16918") == 0 ) {
  //    fprintf(stderr,"Have %s\n",p->name);
  //    int j,seqid;
  //    if ( p->type == BWA_TYPE_NO_MATCH ) {fprintf(stderr," It is unmapped!\n"); j=1;}
  //    else { fprintf(stderr,"  It is mapped.\n");j = pos_end(p) - p->pos;}
  //    bns_coor_pac2real(bns, p->pos, j, &seqid);
  //    if ( p->strand ) fprintf(stderr,"   -");
  //    else fprintf(stderr,"   +");
  //    fprintf(stderr,"%s:%d\n",bns->anns[seqid].name,(int)(p->pos - bns->anns[seqid].offset + 1));
  //  }
   

	int j;
	if (p->type != BWA_TYPE_NO_MATCH || (mate && mate->type != BWA_TYPE_NO_MATCH)) {
		int seqid, nn, am = 0, flag = p->extra_flag;
		char XT;

		if (p->type == BWA_TYPE_NO_MATCH) {
			p->pos = mate->pos;
			p->strand = mate->strand;
			flag |= SAM_FSU;
			j = 1;
		} else j = pos_end(p) - p->pos; // j is the length of the reference in the alignment

		// get seqid
		nn = bns_coor_pac2real(bns, p->pos, j, &seqid);
		if (p->type != BWA_TYPE_NO_MATCH && p->pos + j - bns->anns[seqid].offset > bns->anns[seqid].len) {
			flag |= SAM_FSU; // flag UNMAP as this alignment bridges two adjacent reference sequences
                        //                        if ( strcmp(p->name,"62A73AAXX100903:1:100:10014:16918") == 0 ) {
                        //                          fprintf(stderr,"UN-mapping the read %s : bridges adjacent seqs\n",p->name);
                        //                          fprintf (stderr,"  alignment: ");
                        //                          if ( p->strand ) fprintf(stderr,"-");
                        //                          else fprintf(stderr,"+");
                        //                          fprintf(stderr,"%s:%d\n",bns->anns[seqid].name,(int)(p->pos - bns->anns[seqid].offset + 1));
                        //                        }
                }

		// update flag and print it
		if (p->strand) flag |= SAM_FSR;
		if (mate) {
			if (mate->type != BWA_TYPE_NO_MATCH) {
				if (mate->strand) flag |= SAM_FMR;
			} else flag |= SAM_FMU;
		}
		printf("%s\t%d\t%s\t", p->name, flag, bns->anns[seqid].name);
		printf("%d\t%d\t", (int)(p->pos - bns->anns[seqid].offset + 1), p->mapQ);

		// print CIGAR
		if (p->cigar) {
			for (j = 0; j != p->n_cigar; ++j)
				printf("%d%c", __cigar_len(p->cigar[j]), "MIDS"[__cigar_op(p->cigar[j])]);
		} else if (p->type == BWA_TYPE_NO_MATCH) printf("*");
		else printf("%dM", p->len);

		// print mate coordinate
		if (mate && mate->type != BWA_TYPE_NO_MATCH) {
			int m_seqid, m_is_N;
			long long isize;
			am = mate->seQ < p->seQ? mate->seQ : p->seQ; // smaller single-end mapping quality
			// redundant calculation here, but should not matter too much
			m_is_N = bns_coor_pac2real(bns, mate->pos, mate->len, &m_seqid);
			printf("\t%s\t", (seqid == m_seqid)? "=" : bns->anns[m_seqid].name);
			isize = (seqid == m_seqid)? pos_5(mate) - pos_5(p) : 0;
			if (p->type == BWA_TYPE_NO_MATCH) isize = 0;
			printf("%d\t%lld\t", (int)(mate->pos - bns->anns[m_seqid].offset + 1), isize);
		} else if (mate) printf("\t=\t%d\t0\t", (int)(p->pos - bns->anns[seqid].offset + 1));
		else printf("\t*\t0\t0\t");

		// print sequence and quality
                bwa_print_seq_and_qual(p);

		if (p->clip_len < p->full_len) printf("\tXC:i:%d", p->clip_len);
		if (p->type != BWA_TYPE_NO_MATCH) {
			int i;
			// calculate XT tag
			XT = "NURM"[p->type];
			if (nn > 10) XT = 'N';
			// print tags
			printf("\tXT:A:%c\t%s:i:%d", XT, (mode & BWA_MODE_COMPREAD)? "NM" : "CM", p->nm);
			if (nn) printf("\tXN:i:%d", nn);
			if (mate) printf("\tSM:i:%d\tAM:i:%d", p->seQ, am);
			if (p->type != BWA_TYPE_MATESW) { // X0 and X1 are not available for this type of alignment
				printf("\tX0:i:%d", p->c1);
				if (p->c1 <= max_top2) printf("\tX1:i:%d", p->c2);
			}
			printf("\tXM:i:%d\tXO:i:%d\tXG:i:%d", p->n_mm, p->n_gapo, p->n_gapo+p->n_gape);
			if (p->md) printf("\tMD:Z:%s", p->md);
			// print multiple hits
			if (p->n_multi) {
				printf("\tXA:Z:");
				for (i = 0; i < p->n_multi; ++i) {
					bwt_multi1_t *q = p->multi + i;
					int k;
					j = pos_end_multi(q, p->len) - q->pos;
					nn = bns_coor_pac2real(bns, q->pos, j, &seqid);
					printf("%s,%c%d,", bns->anns[seqid].name, q->strand? '-' : '+',
						   (int)(q->pos - bns->anns[seqid].offset + 1));
					if (q->cigar) {
						for (k = 0; k < q->n_cigar; ++k)
							printf("%d%c", __cigar_len(q->cigar[k]), "MIDS"[__cigar_op(q->cigar[k])]);
					} else printf("%dM", p->len);
					printf(",%d;", q->gap + q->mm);
				}
			}
		}
		putchar('\n');
	} else { // this read has no match
		ubyte_t *s = p->strand? p->rseq : p->seq;
		int flag = p->extra_flag | SAM_FSU;
		if (mate && mate->type == BWA_TYPE_NO_MATCH) flag |= SAM_FMU;
		printf("%s\t%d\t*\t0\t0\t*\t*\t0\t0\t", p->name, flag);
		for (j = 0; j != p->len; ++j) putchar("ACGTN"[(int)s[j]]);
		putchar('\t');
		if (p->qual) {
                  if (p->strand) {
                    //                    seq_reverse(p->len, p->qual, 0); // reverse quality
                    for ( s = p->qual + p->len - 1 ; s >= p->qual ; s-- ) putchar( *s );
                    
                  } else {
			printf("%s", p->qual);
                  }
		} else printf("*");
		if (p->clip_len < p->full_len) printf("\tXC:i:%d", p->clip_len);
		putchar('\n');
	}
}

/* UNUSED
void bwa_print_partial_sam1(const bntseq_t *bns, bwa_seq_t *p, const bwa_seq_t *mate, int mode, int max_top2)
{
	int j;
	if (p->type != BWA_TYPE_NO_MATCH || (mate && mate->type != BWA_TYPE_NO_MATCH)) {
		int seqid, nn, am = 0, flag = p->extra_flag;
		char XT;

		if (p->type == BWA_TYPE_NO_MATCH) {
			p->pos = mate->pos;
			p->strand = mate->strand;
			flag |= SAM_FSU;
			j = 1;
		} else j = pos_end(p) - p->pos; // j is the length of the reference in the alignment

		// get seqid
		nn = bns_coor_pac2real(bns, p->pos, j, &seqid);
		if (p->type != BWA_TYPE_NO_MATCH && p->pos + j - bns->anns[seqid].offset > bns->anns[seqid].len)
			flag |= SAM_FSU; // flag UNMAP as this alignment bridges two adjacent reference sequences

		// update flag and print it
		if (p->strand) flag |= SAM_FSR;
		if (mate) {
			if (mate->type != BWA_TYPE_NO_MATCH) {
				if (mate->strand) flag |= SAM_FMR;
			} else flag |= SAM_FMU;
		}
		printf("%s\t%d\t%s\t", p->name, flag, bns->anns[seqid].name);
		printf("%d\t%d\t", (int)(p->pos - bns->anns[seqid].offset + 1), p->mapQ);

		// print CIGAR
		if (p->cigar) {
			for (j = 0; j != p->n_cigar; ++j)
				printf("%d%c", __cigar_len(p->cigar[j]), "MIDS"[__cigar_op(p->cigar[j])]);
		} else if (p->type == BWA_TYPE_NO_MATCH) printf("*");
		else printf("%dM", p->len);

		// print mate coordinate
		if (mate && mate->type != BWA_TYPE_NO_MATCH) {
			int m_seqid, m_is_N;
			long long isize;
			am = mate->seQ < p->seQ? mate->seQ : p->seQ; // smaller single-end mapping quality
			// redundant calculation here, but should not matter too much
			m_is_N = bns_coor_pac2real(bns, mate->pos, mate->len, &m_seqid);
			printf("\t%s\t", (seqid == m_seqid)? "=" : bns->anns[m_seqid].name);
			isize = (seqid == m_seqid)? pos_5(mate) - pos_5(p) : 0;
			if (p->type == BWA_TYPE_NO_MATCH) isize = 0;
			printf("%d\t%lld\t", (int)(mate->pos - bns->anns[m_seqid].offset + 1), isize);
		} else if (mate) printf("\t=\t%d\t0\t", (int)(p->pos - bns->anns[seqid].offset + 1));
		else printf("\t*\t0\t0\t");


		if (p->clip_len < p->full_len) printf("\tXC:i:%d", p->clip_len);
		if (p->type != BWA_TYPE_NO_MATCH) {
			int i;
			// calculate XT tag
			XT = "NURM"[p->type];
			if (nn > 10) XT = 'N';
			// print tags
			printf("\tXT:A:%c\t%s:i:%d", XT, (mode & BWA_MODE_COMPREAD)? "NM" : "CM", p->nm);
			if (nn) printf("\tXN:i:%d", nn);
			if (mate) printf("\tSM:i:%d\tAM:i:%d", p->seQ, am);
			printf("\tXM:i:%d\tXO:i:%d\tXG:i:%d", p->n_mm, p->n_gapo, p->n_gapo+p->n_gape);
			if (p->md) printf("\tMD:Z:%s", p->md);
		}
		putchar('\n');
	} else { // this read has no match
		ubyte_t *s = p->strand? p->rseq : p->seq;
		int flag = p->extra_flag | SAM_FSU;
		if (mate && mate->type == BWA_TYPE_NO_MATCH) flag |= SAM_FMU;
		printf("%d\t*\t0\t0\t*\t*\t0\t0\t", flag);
		if (p->clip_len < p->full_len) printf("\tXC:i:%d", p->clip_len);
		putchar('\n');
	}
}

*/

bntseq_t *bwa_open_nt(const char *prefix)
{
	bntseq_t *ntbns;
	char *str;
	str = (char*)calloc(strlen(prefix) + 10, 1);
	strcat(strcpy(str, prefix), ".nt");
	ntbns = bns_restore(str);
	free(str);
	return ntbns;
}

void bwa_print_sam_SQ(const bntseq_t *bns)
{
	int i;
	for (i = 0; i < bns->n_seqs; ++i)
		printf("@SQ\tSN:%s\tLN:%d\n", bns->anns[i].name, bns->anns[i].len);
}

void bwase_initialize() 
{
	int i;
	for (i = 1; i != 256; ++i) g_log_n[i] = (int)(4.343 * log(i) + 0.5);
}

void bwa_sai2sam_se_core(const char *prefix, const char *fn_sa, const char *fn_fa, int n_occ)
{
	int i, n_seqs, tot_seqs = 0, m_aln;
	bwt_aln1_t *aln = 0;
	bwa_seq_t *seqs;
	bwa_seqio_t *ks;
	clock_t t;
	bntseq_t *bns, *ntbns = 0;
	FILE *fp_sa;
	gap_opt_t opt;

	// initialization
	bwase_initialize();
	bns = bns_restore(prefix);
	srand48(bns->seed);
	ks = bwa_seq_open(fn_fa);
	fp_sa = xopen(fn_sa, "r");

	// core loop
	m_aln = 0;
	fread(&opt, sizeof(gap_opt_t), 1, fp_sa);
	if (!(opt.mode & BWA_MODE_COMPREAD)) // in color space; initialize ntpac
		ntbns = bwa_open_nt(prefix);
	bwa_print_sam_SQ(bns);
	while ((seqs = bwa_read_seq(ks, 0x40000, &n_seqs, opt.mode & BWA_MODE_COMPREAD, opt.trim_qual)) != 0) {
		tot_seqs += n_seqs;
		t = clock();

		// read alignment
		for (i = 0; i < n_seqs; ++i) {
			bwa_seq_t *p = seqs + i;
			int n_aln;
			fread(&n_aln, 4, 1, fp_sa);
			if (n_aln > m_aln) {
				m_aln = n_aln;
				aln = (bwt_aln1_t*)realloc(aln, sizeof(bwt_aln1_t) * m_aln);
			}
			fread(aln, sizeof(bwt_aln1_t), n_aln, fp_sa);
			bwa_aln2seq_core(n_aln, aln, p, 1, n_occ);
                        //                        seq_reverse(p->len, p->seq, 0);
		}

		fprintf(stderr, "[bwa_aln_core] convert to sequence coordinate... ");
		bwa_cal_pac_pos(prefix, n_seqs, seqs, opt.max_diff, opt.fnr); // forward bwt will be destroyed here
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();

		for (i = 0; i < n_seqs; ++i) {
			bwa_seq_t *p = seqs + i;
                        seq_reverse(p->len, p->seq, 0);
		}

		fprintf(stderr, "[bwa_aln_core] refine gapped alignments... ");
		bwa_refine_gapped(bns, n_seqs, seqs, 0, ntbns);
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();

		fprintf(stderr, "[bwa_aln_core] print alignments... ");
		for (i = 0; i < n_seqs; ++i)
			bwa_print_sam1(bns, seqs + i, 0, opt.mode, opt.max_top2);
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();

		bwa_free_read_seq(n_seqs, seqs);
		fprintf(stderr, "[bwa_aln_core] %d sequences have been processed.\n", tot_seqs);
	}

	// destroy
	bwa_seq_close(ks);
	if (ntbns) bns_destroy(ntbns);
	bns_destroy(bns);
	fclose(fp_sa);
	free(aln);
}

void bwa_print_all_hits(const char *prefix, const char *fn_sa, const char *fn_fa, int max_extra_occ)
{
  int i, n_seqs, tot_seqs = 0, m_aln, m_rest;
	bwt_aln1_t *aln = 0;
	bwa_seq_t *seqs;
	bwa_seqio_t *ks;
	clock_t t,t_convert, t_refine, t_write;;
	bntseq_t *bns, *ntbns = 0;
	FILE *fp_sa;
	gap_opt_t opt;

        //****** below modified (added) for multiple hit printout:

        bwa_seq_t * rest_seqs = 0; // this array will keep (shallow) replicas of the current sequence;
                                   // each of the replicas will be updated with its own alignment
                                   // selected from all the (multiple) alignmens available for the current seq.

        bwt_t *bwt[2];
        char str[1024];
        ubyte_t *pacseq;

	t = clock();
        fprintf(stderr, "[bwa_aln_core] Data structures initialized: ");


	strcpy(str, prefix); strcat(str, ".bwt");  bwt[0] = bwt_restore_bwt(str);
	strcpy(str, prefix); strcat(str, ".sa"); bwt_restore_sa(str, bwt[0]);

	// load reverse BWT and SA
	strcpy(str, prefix); strcat(str, ".rbwt"); bwt[1] = bwt_restore_bwt(str);
	strcpy(str, prefix); strcat(str, ".rsa"); bwt_restore_sa(str, bwt[1]);
        //***************

	// initialization
	bwase_initialize();
	bns = bns_restore(prefix);
	srand48(bns->seed);
	ks = bwa_seq_open(fn_fa);
	fp_sa = xopen(fn_sa, "r");

	pacseq = (ubyte_t*)calloc(bns->l_pac/4+1, 1);
	rewind(bns->fp_pac);
	fread(pacseq, 1, bns->l_pac/4+1, bns->fp_pac);


	// core loop
	m_aln = 0;
        m_rest = 0;

	fread(&opt, sizeof(gap_opt_t), 1, fp_sa);
	if (!(opt.mode & BWA_MODE_COMPREAD)) // in color space; initialize ntpac
		ntbns = bwa_open_nt(prefix);
	bwa_print_sam_SQ(bns);

        fprintf(stderr, "%.2f sec\n", (float)(clock()-t) / CLOCKS_PER_SEC); 

	t = clock();

        max_extra_occ++; // now this variable holds TOTAL number of alignments we want to print (1+requested extra).

	while ((seqs = bwa_read_seq(ks, 0x40000, &n_seqs, opt.mode & BWA_MODE_COMPREAD, opt.trim_qual)) != 0) {
		tot_seqs += n_seqs;
                t_convert = 0;
                t_refine = 0;
                t_write = 0;

                fprintf(stderr, "[bwa_aln_core] %d sequences loaded: ",n_seqs);
                fprintf(stderr, "%.2f sec\n", (float)(clock()-t) / CLOCKS_PER_SEC); 


		// read alignment
		for (i = 0; i < n_seqs; ++i) {

			bwa_seq_t *p = seqs + i;
			int n_aln, n_occ, k, rest;
			fread(&n_aln, 4, 1, fp_sa);
			if (n_aln > m_aln) {
				m_aln = n_aln;
				aln = (bwt_aln1_t*)realloc(aln, sizeof(bwt_aln1_t) * m_aln);
			}

			fread(aln, sizeof(bwt_aln1_t), n_aln, fp_sa);
                        for ( k = n_occ = 0 ; k < n_aln; ++k ) {
				const bwt_aln1_t *q = aln + k;
				n_occ += q->l - q->k + 1;
                        } /* n_occ is now keeping total number of available alignments to the reference 
                             (i.e. placements, NOT bwa records, each of which can describe few placements)
                          */

                        // we are going to keep and print 'rest' alignments:
			rest = ((n_occ > max_extra_occ)? max_extra_occ : n_occ); 

                        if ( rest == 0 ) rest++; /* we need at least one record, even if it is going to say "UNMAPPED" */
                        
                        if ( rest > m_rest ) { 
                          // reallocate rest_seqs array (only if needed) to ensure it can keep 'rest' records
                          m_rest = rest;
                          rest_seqs = (bwa_seq_t*)realloc(rest_seqs,sizeof(bwa_seq_t)*m_rest);
                        }
                        // initialize 'rest' replicas of the current sequence record
                        for ( k = 0 ; k < rest ; k++ ) {
                          rest_seqs[k] = *p; /* clone current sequence p; IMPORTANT: it's a shallow copy */
                        }

                        bwa_aln2seq_all(n_aln, aln, rest,rest_seqs);
                        // now each of the replicas carries its own bwa alignment selected from all alignments
                        // available for the current sequence *p.

                        /* compute positions of the alignments on the ref: */
                        t = clock();

                        bwa_cal_pac_pos_core(bwt[0],bwt[1], rest, rest_seqs, opt.max_diff, opt.fnr ); 
                        t_convert += ( clock() - t );
                        *p = rest_seqs[0]; // copy first selected alignment back into p; 

                        seq_reverse(p->len, p->seq,0);
                
                        /* compute positions of the alignments on the ref: */
                        t = clock();

                        bwa_cal_pac_pos_core(bwt[0],bwt[1], rest, rest_seqs, opt.max_diff, opt.fnr ); 
                        t_convert += ( clock() - t );
                
                        t = clock();

                        bwa_refine_gapped(bns,rest,rest_seqs,pacseq,ntbns); // refine all gapped aligns in our replicas;
                                                                   // side effect: cigars will be allocated for each replica
                        t_refine += ( clock() - t );
                
                        t = clock();
                        //                        for ( k = 0 ; k < n_seqs ; k++ ) {
                        for ( k = 0 ; k < rest ; k++ ) {

                          bwa_print_sam1(bns, rest_seqs + k, 0, opt.mode, opt.max_top2);
                          // cigar was allocated for us in every replica as a side effect, free it now:
                          free ( (rest_seqs+k)->cigar ); 
                        }
                        t_write+= ( clock()-t);
                          
                }

		bwa_free_read_seq(n_seqs, seqs);

		fprintf(stderr, "[bwa_aln_core] convert %d sequences to sequence coordinate: ",n_seqs);
                fprintf(stderr, "%.2f sec\n", (float)t_convert / CLOCKS_PER_SEC); 
		fprintf(stderr, "[bwa_aln_core] refine gapped alignments for %d sequences: ", n_seqs);
		fprintf(stderr, "%.2f sec\n", (float)t_refine / CLOCKS_PER_SEC);
		fprintf(stderr, "[bwa_aln_core] print alignments for %d sequences: ", n_seqs);
		fprintf(stderr, "%.2f sec\n", (float)t_write/ CLOCKS_PER_SEC); 
		fprintf(stderr, "[bwa_aln_core] %d sequences have been processed.\n", tot_seqs);

                t = clock();
	}

	// destroy
	bwt_destroy(bwt[0]);
	bwt_destroy(bwt[1]);

        free(rest_seqs);
        free(pacseq);

	bwa_seq_close(ks);
	if (ntbns) bns_destroy(ntbns);
	bns_destroy(bns);
	fclose(fp_sa);
	free(aln);
}

int bwa_sai2sam_se(int argc, char *argv[])
{
	int c, n_occ = 3;
        int do_full_sam = 0;
	while ((c = getopt(argc, argv, "hsn:f:")) >= 0) {
		switch (c) {
		case 'h': break;
                case 's': do_full_sam = 1; break;
		case 'n': n_occ = atoi(optarg); break;
        case 'f': freopen(optarg, "w", stdout); break;
		default: return 1;
		}
	}

	if (optind + 3 > argc) {
		fprintf(stderr, "Usage: bwa samse [-n max_occ [-s] ] [-f out.sam] <prefix> <in.sai> <in.fq>\n");
		return 1;
	}
        if ( do_full_sam ) bwa_print_all_hits(argv[optind], argv[optind+1], argv[optind+2], n_occ);
	else bwa_sai2sam_se_core(argv[optind], argv[optind+1], argv[optind+2], n_occ);
	return 0;
}
