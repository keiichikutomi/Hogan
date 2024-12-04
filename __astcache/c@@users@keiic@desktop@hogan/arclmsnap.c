/* ========================================================= */
/* PROGRAM GNSHN CR FOR OUTER SURFACE OF LUNAR MARS PROJECT  */
/* DYNAMIC ANALYSIS FOR LARGE DEFORMATION                    */
/* USING CR FORMULATION                                      */
/* CODED BY KEIICHI KUTOMI SINSE 2024.05.26                  */
/* ========================================================= */

int arclmSnap(struct arclmframe* af)
{

	arclmStatic(af);
	arclmDynamic(af);
	return 0;
}

