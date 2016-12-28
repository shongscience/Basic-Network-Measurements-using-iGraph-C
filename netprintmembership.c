/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/ 

#include <igraph.h>
#include <unistd.h>
int main(int argc, char *argv[]) 
{


//	const char *fname = "out.edge";
//	const char *f2name = "out.network";
	igraph_t graph;
	igraph_vector_t weights;

	igraph_vector_t dgresult;
	igraph_vector_t ccresult;
	igraph_vector_t bcresult;
	igraph_vector_t submem;
	igraph_vector_t submemsize;
	igraph_integer_t submemnum;
	igraph_vector_t clresult; //closeness
	igraph_vector_t prresult; //page rank result damp = 0.85
	igraph_vector_t prresultb; //page rank result damp = 0.7
	igraph_vector_t prresultc; //page rank result damp = 0.5

	igraph_vector_t egvec; //eigen vector
	igraph_real_t egval; //eigen value
	igraph_arpack_options_t options;

	igraph_integer_t diaresult;
	igraph_real_t cluresult;
	igraph_real_t diareal;

	int i,ret, vsize,esize;
	FILE *infile;
	FILE *outfile;
	double paircounts;

	// reading arguments
	//printf( "Total input argc = %d\n", argc );
	if ( argc != 3 ) /* argc should be 2 for correct execution */
	{
		//printf( "Argc != 3,  argc = %d\n", argc );
		/* We print argv[0] assuming it is the program name */
		printf( "usage: %s in.edge out.network", argv[0] );
		exit(1);
	}




	/**************************** simple network **************************/
	// read edges and echo them 
	infile=fopen(argv[1], "r");
	if (!infile) {
		printf("Cannot open file: %s\n",argv[1]);
		exit(1);
	}   
	igraph_read_graph_edgelist(&graph, infile,/* max vertices number */ 6 ,/*undirected=*/ 0); 
	ret=fclose(infile);
	if (ret) {
		printf("Cannot close file: %s\n",argv[1]);
		exit(2);
	}   
//	printf("Input Edges :\n");
//	igraph_write_graph_edgelist(&graph, stdout); // echo read-edges
	/////////////////////////////////////////


	//// Diameter
	igraph_diameter(&graph, &diaresult, 0, 0, 0, IGRAPH_UNDIRECTED, 1);
	printf("Diameter: %li\n", (long int) diaresult);
	diareal = (double) diaresult;

	//// Clustering 
	igraph_transitivity_undirected(&graph, &cluresult, IGRAPH_TRANSITIVITY_NAN);
	printf("Clustering: %10f\n", (double) cluresult);

	//// subcomponent ;; plus giant component
	igraph_vector_init(&submem, 0);
	igraph_vector_init(&submemsize, 0);
	igraph_clusters( &graph, &submem, &submemsize,&submemnum,IGRAPH_WEAK);
	printf("Maximum memsize is %10f, vertex %2i.\n",(double)igraph_vector_max(&submemsize), (int)igraph_vector_which_max(&submemsize));
	vsize = (int)igraph_vector_size(&submem); 
	printf("Subcomponent: total node size is %10i\n",vsize);

	




	// write edges and echo them 
	outfile=fopen(argv[2],"w");
    if (!outfile) {
      printf("Cannot open file: %s\n",argv[2]);
      exit(1);
    }
	// write global results with # heading
	fprintf(outfile,"#Diameter: %li\n", (long int) diaresult);
	fprintf(outfile,"#Clustering: %10f\n", (double) cluresult);
	fprintf(outfile,"#GiantComp: %li\n",(long int)igraph_vector_max(&submemsize));
	fprintf(outfile,"#Total node size: %10i\n",vsize);

	// print connected member id
	fprintf(outfile,"#NodeID MemberID\n");
	for (i=0;i < vsize ;i++) {
		fprintf(outfile,"%10i %10i\n",i, (int) VECTOR(submem)[i]);
	}


   
	ret=fclose(outfile);
	if (ret) {
		printf("Cannot close file: %s\n",argv[2]);
		exit(2);
	}   






	// closing
	igraph_vector_destroy(&dgresult);
	igraph_vector_destroy(&ccresult);
	igraph_vector_destroy(&bcresult);
	igraph_vector_destroy(&submem);
	igraph_vector_destroy(&submemsize);
	igraph_vector_destroy(&prresult);
	igraph_vector_destroy(&clresult);
	igraph_vector_destroy(&weights);
	igraph_vector_destroy(&egvec);
//	igraph_integer_destroy(&diaresult);
//	igraph_integer_destory(&cluresult);
	igraph_destroy(&graph);

	return 0;
}
