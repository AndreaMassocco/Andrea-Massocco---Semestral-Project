/************************************************************************************************
*                                                                                               *
*                          Process  module:  Particle Tracking.                                 *
*                                                                                               *
************************************************************************************************/
//#define DEBUG_SW
//#define DEBUG_CONNECTIVITY
//#define DEBUG_INIT
//#define DEBUG_CHARGE
//#define DEBUG_CLUSTER

#include "../hd/po_track.h"

EXTERN int mc_po_track_init()
/*----------------------------------------------------------------------------------------------+
|  This module performs several things:                                                         |
|     1) Define particle tracking plane														    |
|     2) Reference the particles to elements													|
|     3) Track particles                                                                        |
|			3a) 3D Interpolation																|
|			3b) New particle referencing														|
+-----------------------------------------------------------------------------------------------+
|  Return value: 0 if an error occurred.                                                        |
+----------------------------------------------------------------------------------------------*/
{
	double	z = 0, sum_v = 0;
	double	t_while, t_cumulative;
	int		mn = 0, el2f = 0, test = 0, *regions = 0;
	int		*plane_el = 0, el_B = 0, el_t = 0, f_el = 0;
	int		l_plane_el = 0, i = 0, m = 0, j = 0, n = 0, x = 0, l_exit_nodes = 0, yn = 0, s = 0;
	double	area = 0, last_len = -1000;
	int		flag_in = 0, full99 = 0, n_threads = 0;
	int		num_holes = 0;
	FILE   *prof_geom, *debug_test;
	particles *temp_prtcls = 0;

	/*----------------------------define parameters-----------------------------------------------------------*/
	/*--------------------------------------------------------------------------------------------------------*/
	particles_plane = nl_ChargeInterface;	//interface z-position
	n_prtcls = nl_Nprtcl;					//number of particles
	p_qual = nl_pfill;						//the percentual of particles reaching the particles plane
	double	r_corr = 0;					    //correct the radius of the particles interface

	if (nl_SWCWTimeStep > 0)
		t = nl_SWCWTimeStep;
	else
		t = nl_TimeStep;

	save_n = nl_dzsave;					    //distance in mm where the interface is saved

	if (nl_ChargeWeld)
	{
		if (n_prtcls == 0 || p_qual == 0)
		{
			particles_plane = -0.001;
			n_prtcls = 400000;
			save_n = 10;
			printf("\n po: Wrong input values for the Charge Weld read, standard values will be applied.");
			printf("\n po: ChargeInterface = 0");
			printf("\n po: n_prtcls = 400'000");
			printf("\n po: save_n = 10 mm");
		}
	}

	if (nl_SeamWeld)
	{
		save_n = nl_SWdzsave;
		n_prtcls = nl_SWParticlesNumber;
	}

	/*----------------------------allocate variables----------------------------------------------------------*/
	/*--------------------------------------------------------------------------------------------------------*/
	charge_els = (charge_id*)malloc(sizeof(charge_id)*me_nTetras);
	SW_qual_TTrias_K = (double*)malloc(sizeof(double)*rm_nTTrias);
	SW_qual_TTrias_J = (double*)malloc(sizeof(double)*rm_nTTrias);
	id_TTrias_parsing = (int*)malloc(sizeof(int)*rm_nTTrias);
	id_Tnodes_parsing = (int*)malloc(sizeof(int)*rm_nTNodes);
	id_Tetras_parsing = (int*)malloc(sizeof(int)*me_nTetras);
	el_act_regions = (int*)malloc(sizeof(int)*l_plane_el);
	el_regions = (int*)malloc(sizeof(int)*l_plane_el);
	temp_prtcls = (particles*)malloc(sizeof(particles)*n_prtcls);
	start_prtcls = (particles*)malloc(sizeof(particles)*n_prtcls);

	/*--------------------------------------------------------------------------------------------------------*/
	/*Initialize Quality of seam_Welds storage variables*/
	/*--------------------------------------------------------------------------------------------------------*/
	for (i = 0; i < rm_nTTrias; i++)
		id_TTrias_parsing[i] = 0;

	for (i = 0; i < rm_nTNodes; i++)
		id_Tnodes_parsing[i] = 0;

	for (i = 0; i < me_nTetras; i++)
		id_Tetras_parsing[i] = -1;

	/*--------------------------------------------------------------------------------------------------------*/

	/*--------------------------------------------------------------------------------------------------------*/
	/*Initialize Open_MP*/
	/*--------------------------------------------------------------------------------------------------------*/
	if (pa_nProcs_s > 0)
		n_threads = pa_nProcs_s;
	else
		n_threads = omp_get_max_threads();

	omp_set_num_threads(n_threads);

	if (nl_ChargeWeld && save_n <= 0)
	{
		printf("\n po: Negative saving length not allowed. Value corrected to default = 1\n");
		save_n = 1;
	}

	/*--------------------------------------------------------------------------------------------------------*/
	//Find exit_z
	/*--------------------------------------------------------------------------------------------------------*/
	exit_z = -1000;
	for (i = 0; i < co_nTools; i++)
	{
		for (j = 0; j < co_pTool[i].nFacets; j++)
		{
			int m;
			for (m = 0; m < 3; m++)
			{
				int node_cons = co_pTool[i].pFacet[j].Points[m];
				if (co_pTool[i].pPCoord[node_cons][2] > exit_z)
					exit_z = co_pTool[i].pPCoord[node_cons][2];
			}
		}
	}

	if (nl_SeamWeld)
	{
		if (nl_SWQExitProfile > 0)
			exit_z = nl_SWQExitProfile;

		if (nl_SWStartPLane > 0)
			particles_plane = nl_SWStartPLane;
		else
			particles_plane = exit_z / 2;	
	}

	/*--------------------------------------------------------------------------------------------------------*/
	/*Create the particles and distribute them on the charge interface at z*/
	/*--------------------------------------------------------------------------------------------------------*/
	prtcls = get_particles(particles_plane, n_prtcls, r_corr);

	/*--------------------------------------------------------------------------------------------------------*/
	//Isolate the elements that contain the starting particles
	/*--------------------------------------------------------------------------------------------------------*/
	if (nl_ChargeWeld)
	{
		/*Isolate the elements that are near to the plane where the particles lie*/
		if (particles_plane == 0) particles_plane = -0.001;
		plane_el = get_band_el(particles_plane, &l_plane_el);

		el_regions = (int*)realloc(el_regions, sizeof(int)*l_plane_el);
		el_act_regions = (int*)realloc(el_act_regions, sizeof(int)*l_plane_el);

		for (i = 0; i < l_plane_el; i++)
			el_regions[i] = 0;
	}
	else if (nl_SeamWeld)
	{


		plane_el = get_band_el(particles_plane, &l_plane_el);

		el_regions = (int*)realloc(el_regions, sizeof(int)*l_plane_el);
		el_act_regions = (int*)realloc(el_act_regions, sizeof(int)*l_plane_el);

		for (i = 0; i < l_plane_el; i++)
		{
			el_act_regions[i] = 0;
			el_regions[i] = 0;
		}

		connectivity(plane_el, el_act_regions, el_regions, l_plane_el, &num_holes);

		if (num_holes < 2)
		{
			GE_FREE(plane_el);
			int *plane_el;
			particles_plane = exit_z / 4;
			for (i = 0; i < n_prtcls; i++)
				prtcls[i].xyz_new[2] = particles_plane;
			plane_el = get_band_el(particles_plane, &l_plane_el);

			el_act_regions = (int*)realloc(el_act_regions, sizeof(int)*l_plane_el);
			el_regions = (int*)realloc(el_regions, sizeof(int)*l_plane_el);
			for (i = 0; i < l_plane_el; i++)
			{
				el_act_regions[i] = 0;
				el_regions[i] = 0;
			}

			connectivity(plane_el, el_act_regions, el_regions, l_plane_el, &num_holes);
		}

		if (num_holes < 2)
		{
			printf("\n po: No portholes have been found! The last position tried is: %f ", particles_plane);
			printf("\n po: Please check that the inout files are correct or specify a search position through SwStartPLane in the ctl file ");
			printf("\n po: The seam-welds-quality-analysis will be terminated");
			return 0;
		}

		for (i = 0; i < n_prtcls; i++)
		{
			prtcls[i].xyz_new[2] = particles_plane;
			prtcls[i].sw_cons = 0;
		}
	}

	/*--------------------------------------------------------------------------------------------------------*/
	//Initialisation of the particles
	/*--------------------------------------------------------------------------------------------------------*/
	int n_prtcls_ok = 0;
#pragma omp parallel for private(i) reduction(+:n_prtcls_ok)
	for (i = 0; i < n_prtcls; i++)
	{
		check_el(&prtcls[i], plane_el, l_plane_el, el_regions);
		if (prtcls[i].track == 1)
		{
			baricentric(&prtcls[i]);
			avg_v(&prtcls[i]);
			prtcls[i].crossed = 0;
			prtcls[i].qual_sw_K = 0;
			prtcls[i].qual_sw_J = 0;
			prtcls[i].qual_sw_Q = 0;
			prtcls[i].qual_sw_max_p = 0;
			prtcls[i].tracked_length = 0;
			prtcls[i].acc_eps = 0;
			prtcls[i].acc_time = 0;
			n_prtcls_ok++;
		}

	}

#ifdef DEBUG_INIT
	char Name[255];
	debug_test = fopen("charge_out_init.csv", "w");
	fprintf(debug_test, "x,y,z");
	for (j = 0; j < n_prtcls; j++)
		if (prtcls[j].track == 1 && (!isnan(prtcls[j].xyz_new[0]) || !isnan(prtcls[j].xyz_new[1]) || !isnan(prtcls[j].xyz_new[2])))
			fprintf(debug_test, "\n %f,%f,%f", prtcls[j].xyz_new[0], prtcls[j].xyz_new[1], prtcls[j].xyz_new[2]);
	fclose(debug_test);
#endif

	/*--------------------------------------------------------------------------------------------------------*/
	//Delete plane_el elements do not need them anymore
	/*--------------------------------------------------------------------------------------------------------*/
	GE_FREE(plane_el);

	/*--------------------------------------------------------------------------------------------------------*/
	//Take only the particles that are inside the mesh, important for seamwelds
	/*--------------------------------------------------------------------------------------------------------*/
	int k;
	for (k = 0; k < n_prtcls; k++)
		temp_prtcls[k] = prtcls[k];

	prtcls = (particles*)realloc(prtcls, sizeof(particles)*n_prtcls_ok);

	int position = 0;
	for (i = 0; i < n_prtcls; i++)
	{
		if (temp_prtcls[i].track == 1)
		{
			prtcls[position] = temp_prtcls[i];
			position++;
		}
	}
	n_prtcls = n_prtcls_ok;
	GE_FREE(temp_prtcls);

	/*--------------------------------------------------------------------------------------------------------*/
	//Save start situation
	/*--------------------------------------------------------------------------------------------------------*/
	for (i = 0; i < n_prtcls; i++)
		start_prtcls[i] = prtcls[i];

	/*--------------------------------------------------------------------------------------------------------*/
	//Find average velocity at exit
	/*--------------------------------------------------------------------------------------------------------*/
	for (i = 0; i < rm_nTTrias; i++)
	{
		for (j = 0; j < 3; j++)
		{
			l_exit_nodes++;
			sum_v = sum_v + me_pNdata[rm_pTTdata[i].Nodes[j]].Velocity[2];
		}
	}
	exit_v = sum_v / l_exit_nodes;

	/*--------------------------------------------------------------------------------------------------------*/
	//save the exit geometry
	/*--------------------------------------------------------------------------------------------------------*/
	prof_geom = fopen("ProfGeom_Exit.txt", "w");
	for (i = 0; i < rm_nTTrias; i++)
	{
		for (j = 0; j < 3; j++)
		{
			fprintf(prof_geom, "\n %f,%f,%f", me_pNdata[rm_pTTdata[i].Nodes[j]].Coord[0],
				me_pNdata[rm_pTTdata[i].Nodes[j]].Coord[1], me_pNdata[rm_pTTdata[i].Nodes[j]].Coord[2]);
		}
	}
	fclose(prof_geom);

	/*--------------------------------------------------------------------------------------------------------*/
	//Find true outlet area
	/*--------------------------------------------------------------------------------------------------------*/
	A_true = 0;
	for (i = 0; i < rm_nTTrias; i++)
	{
		tria_area(i, &area);
		A_true = A_true + area;
	}

	/*--------------------------------------------------------------------------------------------------------*/
	//Initialise Elements ID
	/*--------------------------------------------------------------------------------------------------------*/
	for (i = 0; i < me_nTetras; i++)
	{
		charge_els[i].element = i;
		charge_els[i].id = 0;
	}

	/*--------------------------------------------------------------------------------------------------------*/
	//Initialize A_mesh
	/*--------------------------------------------------------------------------------------------------------*/
	A_mesh = 0;

	/*--------------------------------------------------------------------------------------------------------*/
	//Find out how many particles are really active
	/*--------------------------------------------------------------------------------------------------------*/
	for (i = 0; i < n_prtcls; i++)
		if (prtcls[i].track == 1)
			n_prtcls_true++;

	next_save_length = 0;

	return 1;
}

EXTERN int mc_po_track()
{
	int search_continue = 1, rep_res = 0, rep_act = 0;
	int n = 0;
	int n_crossed = 0;
	int i = 0, j = 0;
	FILE *charge_out;
	FILE *A_perc;
	int prtcls_crox = 0;
	double n_crossed_last = 0, res_old = 0;
	int no_change = 0;
	/*--------------------------------------------------------------------------------------------------------*/
	//TRACK PARTICLES
	/*--------------------------------------------------------------------------------------------------------*/
	/*--------------------------------------------------------------------------------------------------------*/
	//open text file
	/*--------------------------------------------------------------------------------------------------------*/
	if (nl_ChargeWeld == 1)
	{

		A_perc = fopen("AreaPercentual.txt", "w");
		fprintf(A_perc, "StepNumber,Distance,Area,AreaPercentual,NParticles");
	}

	if (nl_SeamWeld)
		printf("\n po: Welding zones search started");

	
	while (search_continue)
	{
		//PUSH Particles
		track_particle(n_prtcls, n);

		//Find number of crossed particles
		n_crossed = 0; //reset the number of particles crossed
		for (j = 0; j < n_prtcls; j++)
			if (prtcls[j].xyz_new[2] > exit_z)
				n_crossed++;

		if (nl_ChargeWeld)
		{
			//Write to charge_output.csv
			if (rounddbl(exit_v*t*n) == next_save_length)
			{
				//Write in csv file the xyz coordinates of the particles
				char Name[255];

#ifdef DEBUG_CHARGE

				sprintf(Name, "charge_output_debug_%d.csv", rounddbl(exit_v*t*n));
				charge_out = fopen(Name, "w");
				fprintf(charge_out, "x,y,z");
				for (j = 0; j < n_prtcls; j++)
					if (prtcls[j].track == 1 && (!isnan(prtcls[j].xyz_new[0]) || !isnan(prtcls[j].xyz_new[1]) || !isnan(prtcls[j].xyz_new[2])))
						fprintf(charge_out, "\n %f,%f,%f", prtcls[j].xyz_new[0], prtcls[j].xyz_new[1], prtcls[j].xyz_new[2]);
				fclose(charge_out);
#endif
				if (n_crossed > 0)
				{
					sprintf(Name, "charge_output_%d.csv", rounddbl(exit_v*t*n));
					charge_out = fopen(Name, "w");
					fprintf(charge_out, "x,y,z");
					for (j = 0; j < n_prtcls; j++)
						if (prtcls[j].crossed == 1 && prtcls[j].track == 1 && (!isnan(prtcls[j].xyz_new[0]) || !isnan(prtcls[j].xyz_new[1]) || !isnan(prtcls[j].xyz_new[2])))
							fprintf(charge_out, "\n %f,%f,%f", prtcls[j].xyz_new[0], prtcls[j].xyz_new[1], exit_z);
					fclose(charge_out);
				}

				//Check the area covered by new material
				check_trias_area();

				//Write in Area_percentual file
				fprintf(A_perc, "\n %i,%i,%f,%f,%i", n, rounddbl(exit_v*t*n), A_mesh, A_mesh / A_true, n_crossed);

				//Write in output file
				printf("\n po: Charge Weld Analysis progress: %.2f", A_mesh / A_true);

				next_save_length += save_n;

				if (A_mesh / A_true > 0)
				{
					if (A_mesh / A_true == res_old)
					{
						rep_res++;
					}
					else
					{
						res_old = A_mesh / A_true;
					}
				}

				if (rep_res > nl_CWLimitIncrements)
				{
					printf("\n po: Charge Weld Analysis result limit reached. The achieved quality is: %.2f", A_mesh / A_true);
					rep_act = 1;
				}

			}

			// STOP Criterium
			if (A_mesh / A_true >= p_qual || rep_act == 1)
			{
				search_continue = 0;
				write_full_charge(1);
			}

			// CRITICAL STOP Criterium
			if (next_save_length > nl_MaxCWLength)
			{
				printf("\n po: The computation of the charge weld is killed as a MaxCWLength > 50 m achieved ist.");
				printf("\n po: The MaxCWLength can be modified in the ctl File.");
				search_continue = 0;
			}

		}
		else if (nl_SeamWeld)
		{

			//Check which 3D Elements define the welding zone
			for (j = 0; j < n_prtcls; j++)
			{
				int resident_el = prtcls[j].el_id_new;
				if (id_Tetras_parsing[resident_el] == -1)
				{
					id_Tetras_parsing[resident_el] = prtcls[j].id_sw;
				}
				else if (id_Tetras_parsing[resident_el] != -1 && id_Tetras_parsing[resident_el] != prtcls[j].id_sw)
				{
					id_Tetras_parsing[resident_el] = 909;
				}
			}

			if (rounddbl(exit_v*t*n) == next_save_length)
			{
				if (nl_WriteSWTrackToCSV) {
					char Name[255];
					sprintf(Name, "seam_welds_track_%d.csv", rounddbl(exit_v*t*n));
					charge_out = fopen(Name, "w");
					fprintf(charge_out, "id_channel,id,x,y,z,Vx,Vy,Vz,b1,b2,b3,b4,el");
					for (j = 0; j < n_prtcls; j++)
						if (prtcls[j].track == 1 && (!isnan(prtcls[j].xyz_new[0]) && !isnan(prtcls[j].xyz_new[1]) && !isnan(prtcls[j].xyz_new[2])))
						{
							fprintf(charge_out, "\n %i,%i,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%i",
								prtcls[j].id_sw, prtcls[j].id_sw,prtcls[j].xyz_new[0], prtcls[j].xyz_new[1], prtcls[j].xyz_new[2],
								prtcls[j].v_new[0], prtcls[j].v_new[1], prtcls[j].v_new[2],
								prtcls[j].bcoord_new[0], prtcls[j].bcoord_new[1], prtcls[j].bcoord_new[2], prtcls[j].bcoord_new[3],
								prtcls[j].el_id_new);
						}
					fclose(charge_out);
					printf(".");
				}

				if (n_crossed != 0)
					prtcls_crox = 1;

				if (prtcls_crox == 1)
				{
					if (n_crossed != n_crossed_last)
						n_crossed_last = n_crossed;
					else
						no_change++;

					if (no_change > 80)
						search_continue = 0;

					if (n_crossed >= nl_SWpercParticles * n_prtcls_true)
						search_continue = 0;


				}
				next_save_length += save_n;
			}
		}

		//Next integration cycle
		n++;

		//Check how many particles are still moving within the 3D mesh
		int n_prtcls_true = 0;
		for (i = 0; i < n_prtcls; i++)
			if (prtcls[i].track == 1)
				n_prtcls_true++;
	}

	//save the number of increments required to track the particles until nl_SWpercParticles is reached
	number_int_cycles = n;

	if (nl_SeamWeld)
	{
		printf("\n po: Welding zones search completed.");

		//sort 3D elements: give 909 to the elements that have different neighboring indeces 
		for (i = 0; i < me_nTetras; i++)
		{
			int nc = 0;
			if (id_Tetras_parsing[i] > -1 && id_Tetras_parsing[i] < 909)
			{
				int id_sw = id_Tetras_parsing[i];

				int neigh1 = -10000000, neigh2 = -10000000, neigh3 = -100000001, neigh4 = -10000000;

				if (me_pETdata[i].Nbr[0] >= 0)
					neigh1 = id_Tetras_parsing[me_pETdata[i].Nbr[0]];
				if (me_pETdata[i].Nbr[1] >= 0)
					neigh2 = id_Tetras_parsing[me_pETdata[i].Nbr[1]];
				if (me_pETdata[i].Nbr[2] >= 0)
					neigh3 = id_Tetras_parsing[me_pETdata[i].Nbr[2]];
				if (me_pETdata[i].Nbr[3] >= 0)
					neigh4 = id_Tetras_parsing[me_pETdata[i].Nbr[3]];

				if (neigh1 > -1 && neigh1 < 909 && neigh1 != id_sw)
				{
					id_Tetras_parsing[i] = 909;
					id_Tetras_parsing[me_pETdata[i].Nbr[0]] = 909;
				}
				if (neigh2 > -1 && neigh2 < 909 && neigh2 != id_sw)
				{
					id_Tetras_parsing[i] = 909;
					id_Tetras_parsing[me_pETdata[i].Nbr[1]] = 909;
				}
				if (neigh3 > -1 && neigh3 < 909 && neigh3 != id_sw)
				{
					id_Tetras_parsing[i] = 909;
					id_Tetras_parsing[me_pETdata[i].Nbr[2]] = 909;
				}
				if (neigh4 > -1 && neigh4 < 909 && neigh4 != id_sw)
				{
					id_Tetras_parsing[i] = 909;
					id_Tetras_parsing[me_pETdata[i].Nbr[3]] = 909;
				}
			}
		}

		//Assign 909 to the elements when all the elements around him are also 909
		for (i = 0; i < me_nTetras; i++)
		{
			int nc = 0;
			if (id_Tetras_parsing[i] == -1)
			{
				for (j = 0; j < 4; j++)
				{
					if (me_pETdata[i].Nbr[j] >= 0)
						if (id_Tetras_parsing[me_pETdata[i].Nbr[j]] == 909)
							nc++;
				}

				if (nc == 4)
					id_Tetras_parsing[i] = 909;
			}
		}

		//clean up 3D elements standing alone with 909 code
		for (i = 0; i < me_nTetras; i++)
		{
			int nc = 0;
			if (id_Tetras_parsing[i] == 909)
			{
				for (j = 0; j < 4; j++)
				{
					if (me_pETdata[i].Nbr[j] >= 0)
					{
						int idx = 0;
						if (id_Tetras_parsing[me_pETdata[i].Nbr[j]] == 909)
							nc++;
						for (idx = 0; idx < 4; idx++)
						{
							if (me_pETdata[j].Nbr[idx] >= 0)
								if (id_Tetras_parsing[me_pETdata[j].Nbr[idx]] == 909)
									nc++;
						}
					}
				}

				if (nc == 0)
					id_Tetras_parsing[i] = 0;
			}

		}

		for (i = 0; i < me_nTetras; i++)
			charge_els[i].id = id_Tetras_parsing[i];

		if (nl_WriteSWFull)
			write_full_charge(2);

	}
	else if (nl_ChargeWeld)
	{
		fclose(A_perc);
		printf("\n po: Charge weld evolution completed");
	}

	return 1;
}

EXTERN int mc_po_track_sw(int n_iter)
{
	/*--------------------------------------------------------------------------------------------------------*/
	//Extra Routines in case quality of seam welds is required
	/*--------------------------------------------------------------------------------------------------------*/

	int i = 0;
	int j = 0;

	FILE *charge_out;

	/*--------------------------------------------------------------------------------------------------------*/
	//track particles quality
	/*--------------------------------------------------------------------------------------------------------*/
	printf("\n po: Quality analysis of seam welds started ");

	next_save_length = 0;
	for (i = 0; i < me_nTetras; i++)
		charge_els[i].id = id_Tetras_parsing[i];

	for (i = 0; i < n_prtcls; i++)
		prtcls[i] = start_prtcls[i];

	for (j = 0; j < number_int_cycles; j++)
	{
		track_particle(n_prtcls, j);

		for (i = 0; i < n_prtcls; i++)
		{
			int resident_el = prtcls[i].el_id_new;
			double l_track = 0;

			l_track = DIST3(prtcls[i].xyz_new[0], prtcls[i].xyz_new[1], prtcls[i].xyz_new[2],
				prtcls[i].xyz_old[0], prtcls[i].xyz_old[1], prtcls[i].xyz_old[2]);

			if (id_Tetras_parsing[resident_el] == 909 &&
				(prtcls[i].v_new[0] != 0 && prtcls[i].v_new[1] != 0 && prtcls[i].v_new[2] != 0))
			{
				double pressure = 0, yield_stress = 0, eps_rate = 0, temperature = 0;
				int n1 = 0, n2 = 0, n3 = 0, n4 = 0;
				double b1 = 0, b2 = 0, b3 = 0, b4 = 0;
				double R = 8.3144598;
				double Qd = nl_QdSWJ;

				/*if (ma_pMdata[me_pETdata[0].iMd].CurveType == MA_SINEHYPINV)
				Qd = ma_pMdata[0].SineHypInv.QdivbyR * 8.31;
				else if (ma_pMdata[me_pETdata[0].iMd].CurveType == MA_ZENER)
				Qd = ma_pMdata[0].Zener.Q * 8.314;*/

				n1 = me_pETdata[resident_el].Nodes[0];
				n2 = me_pETdata[resident_el].Nodes[1];
				n3 = me_pETdata[resident_el].Nodes[2];
				n4 = me_pETdata[resident_el].Nodes[3];
				b1 = prtcls[i].bcoord_new[0];
				b2 = prtcls[i].bcoord_new[1];
				b3 = prtcls[i].bcoord_new[2];
				b4 = prtcls[i].bcoord_new[3];

				pressure = me_pNdata[n1].Pressure * b1 + me_pNdata[n2].Pressure * b2 + me_pNdata[n3].Pressure * b3 + me_pNdata[n4].Pressure * b4;
				temperature = me_pNdata[n1].Temp * b1 + me_pNdata[n2].Temp * b2 + me_pNdata[n3].Temp * b3 + me_pNdata[n4].Temp * b4;
				eps_rate = me_pETdata[resident_el].dEpqdt;
				yield_stress = fabs(me_pETdata[resident_el].Stressq);

				prtcls[i].tracked_length += l_track;
				prtcls[i].acc_eps += fabs(eps_rate*t);
				prtcls[i].acc_time += t;

				if (pressure < 0)
				{
					pressure = fabs(pressure);
					prtcls[i].qual_sw_K += (pressure / yield_stress)*l_track;
					prtcls[i].qual_sw_J += (pressure / yield_stress)*eps_rate*t*exp(R*temperature / Qd);
					prtcls[i].qual_sw_Q += (pressure / yield_stress)*t;
					if (pressure > prtcls[i].qual_sw_max_p)
						prtcls[i].qual_sw_max_p = pressure;
				}
				else
				{
					if (nl_SWConsiderTensionStresses == 1)
					{
						pressure = fabs(pressure);
						prtcls[i].qual_sw_K -= (pressure / yield_stress)*l_track;
						prtcls[i].qual_sw_J -= (pressure / yield_stress)*eps_rate*t*exp(R*temperature / Qd);
						prtcls[i].qual_sw_Q -= (pressure / yield_stress)*t;
					}
				}
			}
		}

		if (rounddbl(exit_v*t*j) == next_save_length)
		{
			if (nl_WriteSWQualityTrackToCSV) {
				printf("\n po: Seam welds quality progress: %f", (float)i / (float)number_int_cycles * 100);
				char Name[255];
				sprintf(Name, "seam_welds_quality_%d.csv", rounddbl(exit_v*t*j));
				charge_out = fopen(Name, "w");
				fprintf(charge_out, "x,y,z,Vx,Vy,Vz,J,K,Q,max_p");
				int m = 0;
				for (m = 0; m < n_prtcls; m++)
					if (prtcls[m].track == 1 && (!isnan(prtcls[m].xyz_new[0]) && !isnan(prtcls[m].xyz_new[1]) && !isnan(prtcls[m].xyz_new[2])))
					{
						fprintf(charge_out, "\n %f,%f,%f,%f,%f,%f,%f,%f,%f,%f",
							prtcls[m].xyz_new[0], prtcls[m].xyz_new[1], prtcls[m].xyz_new[2],
							prtcls[m].v_new[0], prtcls[m].v_new[1], prtcls[m].v_new[2],
							prtcls[m].qual_sw_J, prtcls[m].qual_sw_K,
							prtcls[m].qual_sw_Q, prtcls[m].qual_sw_max_p);
					}
				fclose(charge_out);
			}
			next_save_length += save_n;
		}
	}

	char Name[255];
	sprintf(Name, "SW_quality_%i.csv", n_iter);
	charge_out = fopen(Name, "w");
	fprintf(charge_out, "ID,x,y,z,J,K,Q,max_p,tracked_length, tracked_time, accumulated_strain");
	int m = 0;
	for (m = 0; m < n_prtcls; m++)
		if (prtcls[m].track == 1 && (!isnan(prtcls[m].xyz_new[0]) && !isnan(prtcls[m].xyz_new[1]) && !isnan(prtcls[m].xyz_new[2])))
		{
			if (prtcls[m].sw_cons == 1 && (prtcls[m].qual_sw_K != 0 || prtcls[m].qual_sw_J != 0 || prtcls[m].qual_sw_Q != 0))
			{
				fprintf(charge_out, "\n %i,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f",
					prtcls[m].id_sw,prtcls[m].xyz_new[0], prtcls[m].xyz_new[1], prtcls[m].xyz_new[2],
					prtcls[m].qual_sw_J, prtcls[m].qual_sw_K,
					prtcls[m].qual_sw_Q, prtcls[m].qual_sw_max_p,
					prtcls[m].tracked_length, prtcls[m].acc_time, prtcls[m].acc_eps);
			}
		}
	fclose(charge_out);
	printf("\n po: Quality analysis of seam welds terminated");

	//Interpolate on nodes
	write_exit_sw_id(n_iter);

	return 1;
}

EXTERN int mc_po_track_release()
{

	GE_FREE(prtcls);
	GE_FREE(charge_els);
	GE_FREE(id_TTrias_parsing);
	GE_FREE(id_Tnodes_parsing);
	GE_FREE(id_Tetras_parsing);
	GE_FREE(start_prtcls);
	GE_FREE(SW_qual_TTrias_K);
	GE_FREE(SW_qual_TTrias_J);
	GE_FREE(el_act_regions);
	GE_FREE(el_regions);
	return 1;
}

/*--------------------------------------------------------------------------------------------------------*/
//SUB-FUNCTIONS
/*--------------------------------------------------------------------------------------------------------*/

void id_charge(charge_id *els, double pplane, int i)
{
	int j, count = 0;
	(*els).element = i;
	for (j = 0; j<4; j++)
	{
		if (me_pNdata[me_pETdata[i].Nodes[j]].Coord[2]<pplane)
		{
			count++;
		}
	}
	if (count == 4)
	{
		(*els).id = 1;
	}
	else if (count == 0)
	{
		(*els).id = 2;
	}
	else
	{
		(*els).id = 0;
	}
}

void avg_v(particles *prtcl)
{

	(*prtcl).v_new[0] = me_pNdata[me_pETdata[(*prtcl).el_id_new].Nodes[0]].Velocity[0] * (*prtcl).bcoord_new[0] +
		me_pNdata[me_pETdata[(*prtcl).el_id_new].Nodes[1]].Velocity[0] * (*prtcl).bcoord_new[1] +
		me_pNdata[me_pETdata[(*prtcl).el_id_new].Nodes[2]].Velocity[0] * (*prtcl).bcoord_new[2] +
		me_pNdata[me_pETdata[(*prtcl).el_id_new].Nodes[3]].Velocity[0] * (*prtcl).bcoord_new[3];

	(*prtcl).v_new[1] = me_pNdata[me_pETdata[(*prtcl).el_id_new].Nodes[0]].Velocity[1] * (*prtcl).bcoord_new[0] +
		me_pNdata[me_pETdata[(*prtcl).el_id_new].Nodes[1]].Velocity[1] * (*prtcl).bcoord_new[1] +
		me_pNdata[me_pETdata[(*prtcl).el_id_new].Nodes[2]].Velocity[1] * (*prtcl).bcoord_new[2] +
		me_pNdata[me_pETdata[(*prtcl).el_id_new].Nodes[3]].Velocity[1] * (*prtcl).bcoord_new[3];

	(*prtcl).v_new[2] = me_pNdata[me_pETdata[(*prtcl).el_id_new].Nodes[0]].Velocity[2] * (*prtcl).bcoord_new[0] +
		me_pNdata[me_pETdata[(*prtcl).el_id_new].Nodes[1]].Velocity[2] * (*prtcl).bcoord_new[1] +
		me_pNdata[me_pETdata[(*prtcl).el_id_new].Nodes[2]].Velocity[2] * (*prtcl).bcoord_new[2] +
		me_pNdata[me_pETdata[(*prtcl).el_id_new].Nodes[3]].Velocity[2] * (*prtcl).bcoord_new[3];
}

void baricentric(particles *prtcl)
{
	(*prtcl).bcoord_new[0] = (*prtcl).d[1] / ((*prtcl).d[0] + 1e-9);
	(*prtcl).bcoord_new[1] = (*prtcl).d[2] / ((*prtcl).d[0] + 1e-9);
	(*prtcl).bcoord_new[2] = (*prtcl).d[3] / ((*prtcl).d[0] + 1e-9);
	(*prtcl).bcoord_new[3] = (*prtcl).d[4] / ((*prtcl).d[0] + 1e-9);
}

void check_el(particles *prtcl, int elems[], int n_els, int el_regions[])
{
	int idx_el, el;

	int n_det;

	double a, b, c, d = 1, e, f, g, h = 1, i, j, k, l = 1, m, n, o, p = 1;

	double det[5];

	for (idx_el = 0; idx_el < n_els; idx_el++)
	{
		el = elems[idx_el];

		for (n_det = 0; n_det < 5; n_det++)
		{

			a = me_pNdata[me_pETdata[el].Nodes[0]].Coord[0];
			b = me_pNdata[me_pETdata[el].Nodes[0]].Coord[1];
			c = me_pNdata[me_pETdata[el].Nodes[0]].Coord[2];

			e = me_pNdata[me_pETdata[el].Nodes[1]].Coord[0];
			f = me_pNdata[me_pETdata[el].Nodes[1]].Coord[1];
			g = me_pNdata[me_pETdata[el].Nodes[1]].Coord[2];

			i = me_pNdata[me_pETdata[el].Nodes[2]].Coord[0];
			j = me_pNdata[me_pETdata[el].Nodes[2]].Coord[1];
			k = me_pNdata[me_pETdata[el].Nodes[2]].Coord[2];

			m = me_pNdata[me_pETdata[el].Nodes[3]].Coord[0];
			n = me_pNdata[me_pETdata[el].Nodes[3]].Coord[1];
			o = me_pNdata[me_pETdata[el].Nodes[3]].Coord[2];


			if (n_det == 1)
			{
				a = prtcl->xyz_new[0];
				b = prtcl->xyz_new[1];
				c = prtcl->xyz_new[2];
			}
			else if (n_det == 2)
			{
				e = prtcl->xyz_new[0];
				f = prtcl->xyz_new[1];
				g = prtcl->xyz_new[2];
			}
			else if (n_det == 3)
			{
				i = prtcl->xyz_new[0];
				j = prtcl->xyz_new[1];
				k = prtcl->xyz_new[2];
			}
			else if (n_det == 4)
			{
				m = prtcl->xyz_new[0];
				n = prtcl->xyz_new[1];
				o = prtcl->xyz_new[2];
			}

			det[n_det] = a*f*k*p - a*f*l*o - a*g*j*p + a*g*l*n + a*h*j*o - a*h*k*n - b*e*k*p + b*e*l*o + b*g*i*p - b*g*l*m - b*h*i*o + b*h*k*m + c*e*j*p - c*e*l*n - c*f*i*p + c*f*l*m + c*h*i*n - c*h*j*m - d*e*j*o + d*e*k*n + d*f*i*o - d*f*k*m - d*g*i*n + d*g*j*m;
		}

		if (fabs(det[0]) < 1E-9) det[0] = 0;
		if (fabs(det[1]) < 1E-9) det[1] = 0;
		if (fabs(det[2]) < 1E-9) det[2] = 0;
		if (fabs(det[3]) < 1E-9) det[3] = 0;
		if (fabs(det[4]) < 1E-9) det[4] = 0;

		if (det[0] != 0 && (det[1] == 0 || det[2] == 0 || det[3] == 0 || det[4] == 0))
		{
			prtcl->d[0] = det[0];
			prtcl->d[1] = det[1];
			prtcl->d[2] = det[2];
			prtcl->d[3] = det[3];
			prtcl->d[4] = det[4];
			prtcl->el_id_new = el;
			break;
		}
		else if (sign(det[0])*sign(det[1]) > 0 && sign(det[0])*sign(det[2]) > 0 &&
			sign(det[0])*sign(det[3]) > 0 && sign(det[0])*sign(det[4]) > 0)
		{
			prtcl->d[0] = det[0];
			prtcl->d[1] = det[1];
			prtcl->d[2] = det[2];
			prtcl->d[3] = det[3];
			prtcl->d[4] = det[4];
			prtcl->el_id_new = el;
			break;
		}
	}

	if (idx_el == n_els)
		prtcl->track = 0;
	else
	{
		prtcl->track = 1;
		prtcl->id_sw = el_regions[idx_el];

	}
	if (prtcl->el_id_new < 0)
		prtcl->track = 0;
}

void check_el2(double x, double y, double z, int elems, int *flag, double *d0, double *d1, double *d2, double *d3, double *d4)
{
	int n_det;

	double a, b, c, d = 1, e, f, g, h = 1, i, j, k, l = 1, m, n, o, p = 1;

	double det[5];

	double p1, p2, p3, p4, p5, p6;

	a = me_pNdata[me_pETdata[elems].Nodes[0]].Coord[0];
	b = me_pNdata[me_pETdata[elems].Nodes[0]].Coord[1];
	c = me_pNdata[me_pETdata[elems].Nodes[0]].Coord[2];

	e = me_pNdata[me_pETdata[elems].Nodes[1]].Coord[0];
	f = me_pNdata[me_pETdata[elems].Nodes[1]].Coord[1];
	g = me_pNdata[me_pETdata[elems].Nodes[1]].Coord[2];

	i = me_pNdata[me_pETdata[elems].Nodes[2]].Coord[0];
	j = me_pNdata[me_pETdata[elems].Nodes[2]].Coord[1];
	k = me_pNdata[me_pETdata[elems].Nodes[2]].Coord[2];

	m = me_pNdata[me_pETdata[elems].Nodes[3]].Coord[0];
	n = me_pNdata[me_pETdata[elems].Nodes[3]].Coord[1];
	o = me_pNdata[me_pETdata[elems].Nodes[3]].Coord[2];

	for (n_det = 0; n_det < 5; n_det++)
	{
		if (n_det == 1)
		{
			a = x;
			b = y;
			c = z;
		}
		else if (n_det == 2)
		{
			a = me_pNdata[me_pETdata[elems].Nodes[0]].Coord[0];
			b = me_pNdata[me_pETdata[elems].Nodes[0]].Coord[1];
			c = me_pNdata[me_pETdata[elems].Nodes[0]].Coord[2];
			e = x;
			f = y;
			g = z;
		}
		else if (n_det == 3)
		{
			e = me_pNdata[me_pETdata[elems].Nodes[1]].Coord[0];
			f = me_pNdata[me_pETdata[elems].Nodes[1]].Coord[1];
			g = me_pNdata[me_pETdata[elems].Nodes[1]].Coord[2];
			i = x;
			j = y;
			k = z;
		}
		else if (n_det == 4)
		{
			i = me_pNdata[me_pETdata[elems].Nodes[2]].Coord[0];
			j = me_pNdata[me_pETdata[elems].Nodes[2]].Coord[1];
			k = me_pNdata[me_pETdata[elems].Nodes[2]].Coord[2];
			m = x;
			n = y;
			o = z;
		}

		/*det[n_det] = a*(f*(k*p - l*o) - g*(j*p - l*n) + h*(j*o - k*n)) -
		b*(e*(k*p - l*o) - g*(i*p - l*m) + h*(i*o - k*m)) +
		c*(e*(j*p - l*n) - f*(i*p - l*m) + h*(i*n - j*m)) -
		d*(e*(j*o - k*n) - f*(i*o - k*m) + g*(i*n - j*m));*/

		p1 = k*p - l*o;
		p2 = j*p - l*n;
		p3 = j*o - k*n;
		p4 = i*p - l*m;
		p5 = i*o - k*m;
		p6 = i*n - j*m;

		det[n_det] = a*(f*p1 - g*p2 + h*p3) -
			b*(e*p1 - g*p4 + h*p5) +
			c*(e*p2 - f*p4 + h*p6) -
			d*(e*p3 - f*p5 + g*p6);

		//det[n_det] = a*f*k*p - a*f*l*o - a*g*j*p + a*g*l*n + a*h*j*o - a*h*k*n - b*e*k*p + b*e*l*o + b*g*i*p - b*g*l*m - b*h*i*o + b*h*k*m + c*e*j*p - c*e*l*n - c*f*i*p + c*f*l*m + c*h*i*n - c*h*j*m - d*e*j*o + d*e*k*n + d*f*i*o - d*f*k*m - d*g*i*n + d*g*j*m;
	}

	if (fabs(det[0]) < 1E-9) det[0] = 0;
	if (fabs(det[1]) < 1E-9) det[1] = 0;
	if (fabs(det[2]) < 1E-9) det[2] = 0;
	if (fabs(det[3]) < 1E-9) det[3] = 0;
	if (fabs(det[4]) < 1E-9) det[4] = 0;

	if (det[0] != 0 && (det[1] == 0 || det[2] == 0 || det[3] == 0 || det[4] == 0))
	{
		*flag = 1;
		*d0 = det[0]; *d1 = det[1]; *d2 = det[2]; *d3 = det[3]; *d4 = det[4];
	}
	else if (sign(det[0])*sign(det[1]) > 0 && sign(det[0])*sign(det[2]) > 0 &&
		sign(det[0])*sign(det[3]) > 0 && sign(det[0])*sign(det[4]) > 0)
	{
		*flag = 1;
		*d0 = det[0]; *d1 = det[1]; *d2 = det[2]; *d3 = det[3]; *d4 = det[4];
	}

}

int sign(double num)
{
	if (num >= 0) return 1;
	else if (num < 0) return -1;
	else return 0;
}

particles* get_particles(double particles_plane, int n_prtcls, double R_corr)
{
	int alpha = 2;
	int i = 0, j;
	double n = n_prtcls;
	double r, theta, x, y, r_tool = 0;
	int b;
	double phi;
	double z = particles_plane;
	double R = me_BilletRadius - me_BilletRadius*0.01 - R_corr;
	particles *prtcls;
	prtcls = (particles*)malloc(n * sizeof(particles));

	/*Find the minimum circle enclosing the die if particles_plane = around 0*/
	if (fabs(particles_plane)<0.1)
	{
		for (i = 0; i < co_nTools; i++)
		{
			for (j = 0; j < co_pTool[i].nPoints; j++)
			{
				if (co_pTool[i].pPCoord[j][2] > 0.1)
					if (pow(pow(co_pTool[i].pPCoord[j][0], 2) + pow(co_pTool[i].pPCoord[j][1], 2), 0.5) > r_tool)
					{
						r_tool = pow(pow(co_pTool[i].pPCoord[j][0], 2) + pow(co_pTool[i].pPCoord[j][1], 2), 0.5);
					}
			}
		}
		r = (R + r_tool) / 2;
		R = r;
	}

	b = (int)(ceil(alpha*sqrt(n)) + 0.1); //to round up 4 sure add 0.1
	phi = (pow(5, 0.5) + 1) / 2;

	for (i = 1; i <= n; i++)
	{
		if (i>n - b) { r = R; }
		else { r = (sqrt(i - 0.5) / (sqrt(n - (b + 1) / 2)))*R; }
		theta = 2 * pi *i / pow(phi, 2);
		x = r*cos(theta);
		y = r*sin(theta);
		prtcls[i - 1].xyz_new[0] = x;
		prtcls[i - 1].xyz_new[1] = y;
		if (me_nSyms>0)
		{
			for (j = 0; j<me_nSyms; j++)
			{
				if (me_pSym[j].Normal[0] != 0)
				{
					prtcls[i - 1].xyz_new[0] = me_pSym[j].Normal[0] * fabs(x);
					if (x == 0) prtcls[i - 1].xyz_new[0] = me_pSym[j].Normal[0] * 0.001;
					continue;
				}
				if (me_pSym[j].Normal[1] != 0)
				{
					prtcls[i - 1].xyz_new[1] = me_pSym[j].Normal[1] * fabs(y);
					if (y == 0) prtcls[i - 1].xyz_new[1] = me_pSym[j].Normal[1] * 0.001;
					continue;
				}
			}
		}
		prtcls[i - 1].xyz_new[2] = z;
		prtcls[i - 1].id = i;
	}
	return prtcls;
}

int *get_band_el(double p_plane, int *l_plane_el)
{
	int i = 0, j = 0, counter = 0, n = 0;
	double n1 = 0, n2 = 0, n3 = 0, n4 = 0;
	int *plane_el;
	double tol = 20;
	for (i = 0; i<me_nElems; i++)
	{
		counter = 0;
		n1 = me_pNdata[me_pETdata[i].Nodes[0]].Coord[2];
		n2 = me_pNdata[me_pETdata[i].Nodes[1]].Coord[2];
		n3 = me_pNdata[me_pETdata[i].Nodes[2]].Coord[2];
		n4 = me_pNdata[me_pETdata[i].Nodes[3]].Coord[2];
		if (n1 <= p_plane) counter++; else counter--;
		if (n2 <= p_plane) counter++; else counter--;
		if (n3 <= p_plane) counter++; else counter--;
		if (n4 <= p_plane) counter++; else counter--;
		if (nl_SeamWeld && abs(counter)<4)
		{
			counter = 5;
			int neigh = -3;
			for (j = 0; j < 4; j++)
			{
				neigh = me_pETdata[i].Nbr[j];
				if (neigh == -2)
				{
					counter = 0;
				}
				else
				{
					int neigh_2_level = -3;
					int m = 0;
					for (m = 0; m < 4; m++)
					{
						neigh_2_level = me_pETdata[neigh].Nbr[m];
						if (neigh_2_level == -2)
						{
							counter = 0;
							break;
						}
						else
						{
							int neigh_3_level = -3;
							int k = 0;
							for (k = 0; k < 4; k++)
							{
								neigh_3_level = me_pETdata[neigh_2_level].Nbr[k];
								if (neigh_3_level == -2)
								{
									counter = 0;
									break;
								}
							}
						}

						if (counter == 0)
							break;
					}
				}

				if (counter == 0)
					break;

			}
		}
		if (abs(counter)<4)
		{
			if (n == 0)
			{
				plane_el = (int*)malloc(sizeof(int));
				plane_el[0] = i;
				n++;
			}
			else
			{
				plane_el = (int*)realloc(plane_el, sizeof(int)*(n + 1));
				plane_el[n] = i;
				*l_plane_el = n + 1;
				n++;
			}
		}
	}
	return plane_el;
}

int get_block_el(particles * prtcl, int *crossed_n1, int *crossed_n2, int *crossed_n3, int *alone_n)
{
	int i = 0, m = 0, flag1 = 0, flag2 = 0, flag3 = 0, found = 0, continue_search = 1, el_idx = -10000000, el_idx_old = -10000000;
	int land_el = -10000000, direction = 1, node1 = -10000000, node2 = -10000000, node3 = -10000000, node_alone = -10000000;
	int el_idx2 = -10000000;
	double delta_z = 0, len_vec = 0;
	double n1x = 0, n1y = 0, n1z = 0, n2x = 0, n2y = 0, n2z = 0, n3x = 0, n3y = 0, n3z = 0, n0x = 0, n0y = 0, n0z = 0;
	double dir[3], P[3], v0[3], v1[3], v2[3];

	double x = prtcl->xyz_new[0], y = prtcl->xyz_new[1], z = prtcl->xyz_new[2];
	double x_old = prtcl->xyz_old[0], y_old = prtcl->xyz_old[1], z_old = prtcl->xyz_old[2];
	double d0, d1, d2, d3, d4;

	el_idx = prtcl->el_id_new;
	el_idx_old = prtcl->el_id_new;
	delta_z = z - z_old;

	/*if (delta_z < 0)
		direction = 0;*/
	direction = 0;

	if (isnan(x))
	{
		int test = 0;
	}

	check_el2(x, y, z, el_idx, &flag1, &d0, &d1, &d2, &d3, &d4);

	if (flag1 == 0)
	{
		while (continue_search)
		{
			flag2 = 0;
			n0x = me_pNdata[me_pETdata[el_idx].Nodes[0]].Coord[0];
			n1x = me_pNdata[me_pETdata[el_idx].Nodes[1]].Coord[0];
			n2x = me_pNdata[me_pETdata[el_idx].Nodes[2]].Coord[0];
			n3x = me_pNdata[me_pETdata[el_idx].Nodes[3]].Coord[0];
			n0y = me_pNdata[me_pETdata[el_idx].Nodes[0]].Coord[1];
			n1y = me_pNdata[me_pETdata[el_idx].Nodes[1]].Coord[1];
			n2y = me_pNdata[me_pETdata[el_idx].Nodes[2]].Coord[1];
			n3y = me_pNdata[me_pETdata[el_idx].Nodes[3]].Coord[1];
			n0z = me_pNdata[me_pETdata[el_idx].Nodes[0]].Coord[2];
			n1z = me_pNdata[me_pETdata[el_idx].Nodes[1]].Coord[2];
			n2z = me_pNdata[me_pETdata[el_idx].Nodes[2]].Coord[2];
			n3z = me_pNdata[me_pETdata[el_idx].Nodes[3]].Coord[2];

			//find normalized direction between starting point and landing point
			len_vec = pow(pow(x - x_old, 2) + pow(y - y_old, 2) + pow(z - z_old, 2), 0.5);

			dir[0] = (x - x_old) / len_vec;
			dir[1] = (y - y_old) / len_vec;
			dir[2] = (z - z_old) / len_vec;

			P[0] = x_old;
			P[1] = y_old;
			P[2] = z_old;

			//find triangles and define verteces	
			for (m = 0; m < 4; m++)
			{
				if (m == 0)
				{
					v0[0] = n1x; v0[1] = n1y; v0[2] = n1z;
					v1[0] = n2x; v1[1] = n2y; v1[2] = n2z;
					v2[0] = n3x; v2[1] = n3y; v2[2] = n3z;
					node1 = me_pETdata[el_idx].Nodes[1];
					node2 = me_pETdata[el_idx].Nodes[2];
					node3 = me_pETdata[el_idx].Nodes[3];
					node_alone = me_pETdata[el_idx].Nodes[0];

				}
				else if (m == 1)
				{
					v0[0] = n0x; v0[1] = n0y; v0[2] = n0z;
					v1[0] = n2x; v1[1] = n2y; v1[2] = n2z;
					v2[0] = n3x; v2[1] = n3y; v2[2] = n3z;
					node1 = me_pETdata[el_idx].Nodes[0];
					node2 = me_pETdata[el_idx].Nodes[2];
					node3 = me_pETdata[el_idx].Nodes[3];
					node_alone = me_pETdata[el_idx].Nodes[1];
				}
				else if (m == 2)
				{
					v0[0] = n0x; v0[1] = n0y; v0[2] = n0z;
					v1[0] = n1x; v1[1] = n1y; v1[2] = n1z;
					v2[0] = n3x; v2[1] = n3y; v2[2] = n3z;
					node1 = me_pETdata[el_idx].Nodes[0];
					node2 = me_pETdata[el_idx].Nodes[1];
					node3 = me_pETdata[el_idx].Nodes[3];
					node_alone = me_pETdata[el_idx].Nodes[2];
				}
				else if (m == 3)
				{
					v0[0] = n0x; v0[1] = n0y; v0[2] = n0z;
					v1[0] = n1x; v1[1] = n1y; v1[2] = n1z;
					v2[0] = n2x; v2[1] = n2y; v2[2] = n2z;
					node1 = me_pETdata[el_idx].Nodes[0];
					node2 = me_pETdata[el_idx].Nodes[1];
					node3 = me_pETdata[el_idx].Nodes[2];
					node_alone = me_pETdata[el_idx].Nodes[3];
				}

				//ray through triangle
				flag2 = rayIntersectsTriangle(P, dir, v0, v1, v2, &direction, len_vec);

				if (direction == 1 && flag2 == 1)
				{
					el_idx2 = me_pETdata[el_idx].Nbr[m];
					if (el_idx2>0)
					{
						check_el2(x, y, z, el_idx2, &flag3, &d0, &d1, &d2, &d3, &d4);
						if (flag3 == 0)
						{
							flag2 = 0;
						}
					}
				}

				if (flag2 == 1)
				{
					if (me_pETdata[el_idx].Nbr[m] != el_idx_old)
					{
						*crossed_n1 = node1;
						*crossed_n2 = node2;
						*crossed_n3 = node3;
						*alone_n = node_alone;
						el_idx_old = el_idx;
						break;
					}
				}
			}

			if (flag2 == 1)
			{
				if (me_pETdata[el_idx].Nbr[m] > 0)
				{
					el_idx = me_pETdata[el_idx].Nbr[m];
					check_el2(x, y, z, el_idx, &flag3, &d0, &d1, &d2, &d3, &d4);

					if (flag3 == 1)
					{
						found = 1;
						continue_search = 0;
					}
				}
				else
				{
					found = 0;
					continue_search = 0;
				}
			}
			else
			{
				found = 0;
				continue_search = 0;
			}

			i++;

			if (i == 500)
			{
				found = 0;
				continue_search = 0;
				el_idx = -el_idx_old;
			}
		}
	}
	else
	{
		found = 1;
	}

	if (found == 1)
	{
		land_el = el_idx;
		prtcl->d[0] = d0;
		prtcl->d[1] = d1;
		prtcl->d[2] = d2;
		prtcl->d[3] = d3;
		prtcl->d[4] = d4;
	}
	else
	{
		land_el = -el_idx;
	}

	return land_el;
}

double max_el_length()
{
	int i, j = 0;
	double l1, l2, l3, l4, l5, l6, l_max;
	double l_max_tot = 0;
	int n1, n2, n3, n4;
	for (i = 0; i < me_nElems; i++)
	{
		n1 = me_pETdata[i].Nodes[0];
		n2 = me_pETdata[i].Nodes[1];
		n3 = me_pETdata[i].Nodes[2];
		n4 = me_pETdata[i].Nodes[3];
		l1 = ndists(n1, n2);
		l2 = ndists(n1, n3);
		l3 = ndists(n1, n4);
		l4 = ndists(n2, n3);
		l5 = ndists(n3, n4);
		l6 = ndists(n2, n4);
		l_max = l1;
		if (l2>l_max) { l_max = l2; }
		if (l3>l_max) { l_max = l3; }
		if (l4>l_max) { l_max = l4; }
		if (l5>l_max) { l_max = l5; }
		if (l6>l_max) { l_max = l6; }
		if (l_max > l_max_tot) { l_max_tot = l_max; }
	}
	return(l_max_tot);
}

double ndists(int n1, int n2)
{
	double x1, x2, y1, y2, z1, z2;
	double dist;
	x1 = me_pNdata[n1].Coord[0];
	x2 = me_pNdata[n2].Coord[0];
	y1 = me_pNdata[n1].Coord[1];
	y2 = me_pNdata[n2].Coord[1];
	z1 = me_pNdata[n1].Coord[2];
	z2 = me_pNdata[n2].Coord[2];
	dist = sqrt(pow((x2 - x1), 2) + pow((y2 - y1), 2) + pow((z2 - z1), 2));
	return(dist);
}

int rayIntersectsTriangle(double *p, double *d,
	double *v0, double *v1, double *v2, int *dir, double length) {

	double e1[3], e2[3], h[3], s[3], q[3];
	double a, f, u, v, t;

	vector(e1, v1, v0);
	vector(e2, v2, v0);

	crossProduct(h, d, e2);
	a = innerProduct(e1, h);

	if (a > -0.000001 && a < 0.000001)
		return 0;

	f = 1 / a;
	vector(s, p, v0);
	u = f * (innerProduct(s, h));

	if (u < 0.0 || u > 1.0)
		return 0;

	crossProduct(q, s, e1);
	v = f * innerProduct(d, q);

	if (v < 0.0 || u + v > 1.0)
		return 0;

	// at this stage we can compute t to find out where
	// the intersection point is on the line
	t = f * innerProduct(e2, q);

	*dir = 0;

	if (fabs(t) < 0.000001)
	{
		t = 0;
		*dir = 1;
	}

	if (t > fabs(length))
		return 0;

	if (t < 0)
		return 0;

	/*if (dir == 1)
	{
		if (t<0)
			return 0;
	}
	if (dir == 0)
	{
		if (t >= 0)
			return 0;
	}*/

	return 1;
}

void baricentric2D(int el, double px, double py, int*yn)
{
	double n1x, n1y, n2x, n2y, n3x, n3y, alpha, beta, gamma;
	n1x = me_pNdata[rm_pTTdata[el].Nodes[0]].Coord[0];
	n1y = me_pNdata[rm_pTTdata[el].Nodes[0]].Coord[1];
	n2x = me_pNdata[rm_pTTdata[el].Nodes[1]].Coord[0];
	n2y = me_pNdata[rm_pTTdata[el].Nodes[1]].Coord[1];
	n3x = me_pNdata[rm_pTTdata[el].Nodes[2]].Coord[0];
	n3y = me_pNdata[rm_pTTdata[el].Nodes[2]].Coord[1];
	*yn = 0;
	alpha = ((n2y - n3y)*(px - n3x) + (n3x - n2x)*(py - n3y)) /
		((n2y - n3y)*(n1x - n3x) + (n3x - n2x)*(n1y - n3y));
	beta = ((n3y - n1y)*(px - n3x) + (n1x - n3x)*(py - n3y)) /
		((n2y - n3y)*(n1x - n3x) + (n3x - n2x)*(n1y - n3y));
	gamma = 1 - alpha - beta;
	if (alpha > 0 && beta>0 && gamma>0)
		*yn = 1;
}

void tria_area(int el, double*area)
{
	double n1x, n1y, n2x, n2y, n3x, n3y;
	n1x = me_pNdata[rm_pTTdata[el].Nodes[0]].Coord[0];
	n1y = me_pNdata[rm_pTTdata[el].Nodes[0]].Coord[1];
	n2x = me_pNdata[rm_pTTdata[el].Nodes[1]].Coord[0];
	n2y = me_pNdata[rm_pTTdata[el].Nodes[1]].Coord[1];
	n3x = me_pNdata[rm_pTTdata[el].Nodes[2]].Coord[0];
	n3y = me_pNdata[rm_pTTdata[el].Nodes[2]].Coord[1];
	*area = fabs(0.5*(n1x*(n2y - n3y) + n2x*(n3y - n1y) + n3x*(n1y - n2y)));
}

int rounddbl(double x)
{
	if (x < 0.0)
		return (int)(x - 0.5);
	else
		return (int)(x + 0.5);
}

void normal2tria(int el, int n_prtcl, double vx, double vy, double vz,
	double *nor_face_x, double *nor_face_y, double *nor_face_z)
{
	double n1x, n1y, n1z, n2x, n2y, n2z, n3x, n3y, n3z;
	double PQ[3], PR[3], normal[3], len, normal2[3];
	double v0f[3], v1f[3], v2f[3], pf[3], dc1, dc2;

	n1x = me_pNdata[prtcls[n_prtcl].tria_crossed[0]].Coord[0];
	n2x = me_pNdata[prtcls[n_prtcl].tria_crossed[1]].Coord[0];
	n3x = me_pNdata[prtcls[n_prtcl].tria_crossed[2]].Coord[0];
	n1y = me_pNdata[prtcls[n_prtcl].tria_crossed[0]].Coord[1];
	n2y = me_pNdata[prtcls[n_prtcl].tria_crossed[1]].Coord[1];
	n3y = me_pNdata[prtcls[n_prtcl].tria_crossed[2]].Coord[1];
	n1z = me_pNdata[prtcls[n_prtcl].tria_crossed[0]].Coord[2];
	n2z = me_pNdata[prtcls[n_prtcl].tria_crossed[1]].Coord[2];
	n3z = me_pNdata[prtcls[n_prtcl].tria_crossed[2]].Coord[2];

	/*ray through triangle*/
	//find direction between starting point and landing point

	v0f[0] = n1x; v0f[1] = n1y; v0f[2] = n1z;
	v1f[0] = n2x; v1f[1] = n2y; v1f[2] = n2z;
	v2f[0] = n3x; v2f[1] = n3y; v2f[2] = n3z;

	pf[0] = me_pNdata[prtcls[n_prtcl].node_alone].Coord[0];
	pf[1] = me_pNdata[prtcls[n_prtcl].node_alone].Coord[1];
	pf[2] = me_pNdata[prtcls[n_prtcl].node_alone].Coord[2];

	PQ[0] = v0f[0] - v1f[0];
	PQ[1] = v0f[1] - v1f[1];
	PQ[2] = v0f[2] - v1f[2];
	PR[0] = v0f[0] - v2f[0];
	PR[1] = v0f[1] - v2f[1];
	PR[2] = v0f[2] - v2f[2];

	normal[0] = PQ[1] * PR[2] - PQ[2] * PR[1];
	normal[1] = PQ[2] * PR[0] - PQ[0] * PR[2];
	normal[2] = PQ[0] * PR[1] - PQ[1] * PR[0];

	normal2[0] = -PQ[1] * PR[2] + PQ[2] * PR[1];
	normal2[1] = -PQ[2] * PR[0] + PQ[0] * PR[2];
	normal2[2] = -PQ[0] * PR[1] + PQ[1] * PR[0];

	len = pow(pow(normal[0], 2) + pow(normal[1], 2) + pow(normal[2], 2), 0.5);

	normal[0] = normal[0] / len;
	normal[1] = normal[1] / len;
	normal[2] = normal[2] / len;
	normal2[0] = normal2[0] / len;
	normal2[1] = normal2[1] / len;
	normal2[2] = normal2[2] / len;

	dc1 = pow(pow(v0f[0] + normal[0] - pf[0], 2) + pow(v0f[1] + normal[1] - pf[1], 2) + pow(v0f[2] + normal[2] - pf[2], 2), 0.5);
	dc2 = pow(pow(v0f[0] + normal2[0] - pf[0], 2) + pow(v0f[1] + normal2[1] - pf[1], 2) + pow(v0f[2] + normal2[2] - pf[2], 2), 0.5);

	if (dc1<dc2)
	{
		*nor_face_x = normal[0];
		*nor_face_y = normal[1];
		*nor_face_z = normal[2];
	}
	else
	{
		*nor_face_x = normal2[0];
		*nor_face_y = normal2[1];
		*nor_face_z = normal2[2];
	}
}

void proj_vec(double *nor_face_x, double *nor_face_y, double *nor_face_z,
	double *proj_v_x, double *proj_v_y, double *proj_v_z,
	double vx, double vy, double vz)
{
	double len, v1, v2, v3, b1, b2, b3;
	len = pow(pow(vx, 2) + pow(vy, 2) + pow(vz, 2), 0.5);
	len = 1;
	v1 = vx / len;
	v2 = vy / len;
	v3 = vz / len;
	b1 = *nor_face_x;
	b2 = *nor_face_y;
	b3 = *nor_face_z;

	double v[3], n[3], c[3], d[3], mag_c = 1;

	v[0] = v1;
	v[1] = v2;
	v[2] = v3;

	n[0] = b1;
	n[1] = b2;
	n[2] = b3;

	crossProduct(c, v, n);
	//mag_c = SQR2(PWR2(c[0]) + PWR2(c[1]) + PWR2(c[2]));

	c[0] = c[0] / mag_c;
	c[1] = c[1] / mag_c;
	c[2] = c[2] / mag_c;

	crossProduct(d, n, c);

	*proj_v_x = d[0] / mag_c;
	*proj_v_y = d[1] / mag_c;
	*proj_v_z = d[2] / mag_c;



	/**proj_v_x = (v1 - ((v1*b1 + v2*b2 + v3*b3) / (b1*b1 + b2*b2 + b3*b3))*b1)*len;
	*proj_v_y = (v2 - ((v1*b1 + v2*b2 + v3*b3) / (b1*b1 + b2*b2 + b3*b3))*b2)*len;
	*proj_v_z = (v3 - ((v1*b1 + v2*b2 + v3*b3) / (b1*b1 + b2*b2 + b3*b3))*b3)*len;*/
	/**proj_v_x = (v1 - (v1*b1 + v2*b2 + v3*b3)*b1)  *len;
	*proj_v_y = (v2 - (v1*b1 + v2*b2 + v3*b3)*b2)  *len;
	*proj_v_z = (v3 - (v1*b1 + v2*b2 + v3*b3)*b3)  *len;*/
}

int rayIntersectsTriangle_back(double *p, double *d,
	double *v0, double *v1, double *v2, int dir, double *t) {

	double e1[3], e2[3], h[3], s[3], q[3];
	double a, f, u, v;
	vector(e1, v1, v0);
	vector(e2, v2, v0);

	crossProduct(h, d, e2);
	a = innerProduct(e1, h);

	if (a > -0.000001 && a < 0.000001)
		return 0;

	f = 1 / a;
	vector(s, p, v0);
	u = f * (innerProduct(s, h));

	if (u < 0.0 || u > 1.0)
		return 0;

	crossProduct(q, s, e1);
	v = f * innerProduct(d, q);

	if (v < 0.0 || u + v > 1.0)
		return 0;

	// at this stage we can compute t to find out where
	// the intersection point is on the line
	*t = f * innerProduct(e2, q);
	if (dir == 1)
	{
		if (*t<0)
			return 0;
	}
	if (dir == 0)
	{
		if (*t >= 0)
			return 0;
	}


	return 1;
}

void backintria(particles *prtcl)
{
	double n1x, n1y, n1z, n2x, n2y, n2z, n3x, n3y, n3z, n4x, n4y, n4z, t1;
	double len_vec, x, x_old, y, y_old, z, z_old, dir[3], P[3], v0[3], v1[3], v2[3], t;
	double v0f[3], v1f[3], v2f[3];
	int m, flag3, el;
	el = (*prtcl).el_id_new;
	n1z = me_pNdata[me_pETdata[el].Nodes[0]].Coord[2];
	n2z = me_pNdata[me_pETdata[el].Nodes[1]].Coord[2];
	n3z = me_pNdata[me_pETdata[el].Nodes[2]].Coord[2];
	n4z = me_pNdata[me_pETdata[el].Nodes[3]].Coord[2];
	n1x = me_pNdata[me_pETdata[el].Nodes[0]].Coord[0];
	n2x = me_pNdata[me_pETdata[el].Nodes[1]].Coord[0];
	n3x = me_pNdata[me_pETdata[el].Nodes[2]].Coord[0];
	n4x = me_pNdata[me_pETdata[el].Nodes[3]].Coord[0];
	n1y = me_pNdata[me_pETdata[el].Nodes[0]].Coord[1];
	n2y = me_pNdata[me_pETdata[el].Nodes[1]].Coord[1];
	n3y = me_pNdata[me_pETdata[el].Nodes[2]].Coord[1];
	n4y = me_pNdata[me_pETdata[el].Nodes[3]].Coord[1];
	x = (*prtcl).xyz_new[0];
	y = (*prtcl).xyz_new[1];
	z = (*prtcl).xyz_new[2];
	x_old = (*prtcl).xyz_old[0];
	y_old = (*prtcl).xyz_old[1];
	z_old = (*prtcl).xyz_old[2];
	/*ray through triangle*/
	//find direction between starting point and landing point
	len_vec = pow(pow(x - x_old, 2) + pow(y - y_old, 2) + pow(z - z_old, 2), 0.5);
	dir[0] = (x - x_old) / len_vec; dir[1] = (y - y_old) / len_vec; dir[2] = (z - z_old) / len_vec;
	P[0] = x_old; P[1] = y_old; P[2] = z_old;
	//find triangles and define verteces	
	for (m = 0; m<4; m++)
	{
		if (m == 0)
		{
			v0[0] = n1x; v0[1] = n1y; v0[2] = n1z;
			v1[0] = n2x; v1[1] = n2y; v1[2] = n2z;
			v2[0] = n3x; v2[1] = n3y; v2[2] = n3z;
			flag3 = rayIntersectsTriangle_back(P, dir, v0, v1, v2, 1, &t);
			t1 = t;
			v0f[0] = n1x; v0f[1] = n1y; v0f[2] = n1z;
			v1f[0] = n2x; v1f[1] = n2y; v1f[2] = n2z;
			v2f[0] = n3x; v2f[1] = n3y; v2f[2] = n3z;
		}
		else if (m == 1)
		{
			v0[0] = n1x; v0[1] = n1y; v0[2] = n1z;
			v1[0] = n3x; v1[1] = n3y; v1[2] = n3z;
			v2[0] = n4x; v2[1] = n4y; v2[2] = n4z;
			flag3 = rayIntersectsTriangle_back(P, dir, v0, v1, v2, 1, &t);
			if (fabs(t)<fabs(t1))
			{
				t1 = t;
				v0f[0] = n1x; v0f[1] = n1y; v0f[2] = n1z;
				v1f[0] = n3x; v1f[1] = n3y; v1f[2] = n3z;
				v2f[0] = n4x; v2f[1] = n4y; v2f[2] = n4z;
			}
		}
		else if (m == 2)
		{
			v0[0] = n1x; v0[1] = n1y; v0[2] = n1z;
			v1[0] = n2x; v1[1] = n2y; v1[2] = n2z;
			v2[0] = n4x; v2[1] = n4y; v2[2] = n4z;
			flag3 = rayIntersectsTriangle_back(P, dir, v0, v1, v2, 1, &t);
			if (fabs(t)<fabs(t1))
			{
				t1 = t;
				v0f[0] = n1x; v0f[1] = n1y; v0f[2] = n1z;
				v1f[0] = n2x; v1f[1] = n2y; v1f[2] = n2z;
				v2f[0] = n4x; v2f[1] = n4y; v2f[2] = n4z;
			}
		}
		else if (m == 3)
		{
			v0[0] = n2x; v0[1] = n2y; v0[2] = n2z;
			v1[0] = n3x; v1[1] = n3y; v1[2] = n3z;
			v2[0] = n4x; v2[1] = n4y; v2[2] = n4z;
			flag3 = rayIntersectsTriangle_back(P, dir, v0, v1, v2, 1, &t);
			if (fabs(t)<fabs(t1))
			{
				t1 = t;
				v0f[0] = n2x; v0f[1] = n1y; v0f[2] = n1z;
				v1f[0] = n3x; v1f[1] = n2y; v1f[2] = n2z;
				v2f[0] = n4x; v2f[1] = n4y; v2f[2] = n4z;
			}
		}
	}
	if (t1<0)
	{
		(*prtcl).xyz_old[0] = x_old + (t1 - 0.01)*dir[0];
		(*prtcl).xyz_old[1] = y_old + (t1 - 0.01)*dir[1];
		(*prtcl).xyz_old[2] = z_old + (t1 - 0.01)*dir[2];
	}
}

EXTERN void track_particle(int n_prtcls, int n)
{
	int		el_B = 0, el_t = 0, i = 0, yn = 0, flag_in = 0, j = 0;
	double	area = 0, A_true = 0;
	double  nor_face_x = 0, nor_face_y = 0, nor_face_z = 0, proj_v_x = 0, proj_v_y = 0, proj_v_z = 0;
	double	temp_x = 0, temp_y = 0, temp_z = 0;


#pragma omp parallel 
#pragma omp for private(j,i,el_B,el_t,temp_x,temp_y,temp_z,flag_in,nor_face_x,nor_face_y,nor_face_z,proj_v_x, proj_v_y, proj_v_z,prtcl_B,yn) schedule(static)
	for (j = 0; j < n_prtcls; j++)
	{
		int crossed_n1 = -10000000, crossed_n2 = -10000000, crossed_n3 = -10000000, alone_n = -10000000;

		if (prtcls[j].track == 1)
		{
			if (fabs(prtcls[j].v_new[0]) < 0.00000001 &&
				fabs(prtcls[j].v_new[1]) < 0.00000001 &&
				fabs(prtcls[j].v_new[2]) < 0.00000001 &&
				prtcls[j].xyz_new[2]<exit_z)
			{
				prtcls[j].track = 0;
			}
			//check if particle is strongly out of geometry
			else if (fabs(prtcls[j].xyz_new[0]) > me_BilletRadius*1.25 || fabs(prtcls[j].xyz_new[1]) > me_BilletRadius*1.25)
			{
				prtcls[j].v_new[0] = 0;
				prtcls[j].v_new[1] = 0;
				prtcls[j].v_new[2] = 0;
				prtcls[j].track = 0;
			}
			//or the particle already left the euler mesh
			else if (prtcls[j].xyz_old[2] + prtcls[j].v_new[2] * t >= exit_z || prtcls[j].xyz_new[2] >= exit_z)
			{
				if (nl_ChargeWeld)
				{
					//set vx vy to 0 while continue vz 
					prtcls[j].v_new[0] = 0;
					prtcls[j].v_new[1] = 0;
					prtcls[j].v_new[2] = exit_v;

					//copy particles new_coord into old_coord
					for (i = 0; i < 3; i++)  prtcls[j].xyz_old[i] = prtcls[j].xyz_new[i];

					//calculate new coordinates in z direction
					prtcls[j].xyz_new[2] = prtcls[j].xyz_old[2] + exit_v*t;
				}
				else if (nl_SeamWeld)
				{
					//set vx vy to 0 while continue vz 
					prtcls[j].v_new[0] = 0;
					prtcls[j].v_new[1] = 0;
					prtcls[j].v_new[2] = 0;
				}

				//check which element is the last crossed element
				if (rounddbl(exit_v*t*n) == next_save_length)
				{
					if (prtcls[j].crossed == 0)
					{
						for (i = 0; i < rm_nTTrias; i++)
						{
							baricentric2D(i, prtcls[j].xyz_new[0], prtcls[j].xyz_new[1], &yn);
							if (yn == 1)
							{
								rm_pTTdata[i].Area_cons = 1;

								if (nl_SeamWeld == 1 && prtcls[j].sw_cons == 0)
								{
									prtcls[j].sw_cons = 1;
									if (id_TTrias_parsing[i] == 0)
										id_TTrias_parsing[i] = prtcls[j].id_sw;
									else if (id_TTrias_parsing[i] != prtcls[j].id_sw && id_TTrias_parsing[i] % 1000 != prtcls[j].id_sw)
										if (id_TTrias_parsing[i] / 1000 != prtcls[j].id_sw)
											id_TTrias_parsing[i] = id_TTrias_parsing[i] * 1000 + prtcls[j].id_sw;

									prtcls[j].el_sw = i;
								}
								break;
							}
						}
						if (nl_ChargeWeld)
							prtcls[j].crossed = 1;
					}
				}
			}
			//or the particle is still inside the euler mesh
			else
			{
				/*check Velocity 2 at half Delta t = Position B*/
				// copy the old coordinates into the B structure
				for (i = 0; i < 3; i++) prtcl_B.xyz_old[i] = prtcls[j].xyz_new[i];

				//calculate coordinates xyz_new for B
				for (i = 0; i < 3; i++) prtcl_B.xyz_new[i] = prtcl_B.xyz_old[i] + prtcls[j].v_new[i] * (t / 2);

				prtcl_B.el_id_new = prtcls[j].el_id_new;
				prtcl_B.el_id_old = prtcls[j].el_id_old;

				//find in which element it lies
				el_B = get_block_el(&prtcl_B, &crossed_n1, &crossed_n2, &crossed_n3, &alone_n);

				//copy particles new velocity into old_velocity
				for (i = 0; i < 3; i++) prtcls[j].v_old[i] = prtcls[j].v_new[i];

				if (el_B >= 0)
				{
					//assign new element
					prtcl_B.el_id_new = el_B;

					//calculate the barycentric coordinates for particle B @ xyz_new
					baricentric(&prtcl_B);

					//calculate new velocity at position B
					avg_v(&prtcl_B);

					//copy new velocity from B
					for (i = 0; i < 3; i++) prtcls[j].v_new[i] = prtcl_B.v_new[i];
				}

				//copy particles new_coord into old_coord
				for (i = 0; i < 3; i++) prtcls[j].xyz_old[i] = prtcls[j].xyz_new[i];

				//calculate new coordinates
				for (i = 0; i < 3; i++) prtcls[j].xyz_new[i] = prtcls[j].xyz_old[i] + prtcls[j].v_new[i] * t;

				/*crossed_n1 = prtcls[j].tria_crossed[0];
				crossed_n2 = prtcls[j].tria_crossed[1];
				crossed_n3 = prtcls[j].tria_crossed[2];
				alone_n = prtcls[j].node_alone;*/

				//find new element
				el_t = get_block_el(&prtcls[j], &crossed_n1, &crossed_n2, &crossed_n3, &alone_n);

				prtcls[j].tria_crossed[0] = crossed_n1;
				prtcls[j].tria_crossed[1] = crossed_n2;
				prtcls[j].tria_crossed[2] = crossed_n3;
				prtcls[j].node_alone = alone_n;

				/*Velocity exclusion for outer boundaries elements*/
				if (el_t < 0)
				{
					double d0, d1, d2, d3, d4;

					temp_x = prtcls[j].xyz_old[0];
					temp_y = prtcls[j].xyz_old[1];
					temp_z = prtcls[j].xyz_old[2];

					//check if last increment was in the element
					check_el2(prtcls[j].xyz_old[0], prtcls[j].xyz_old[1], prtcls[j].xyz_old[2], prtcls[j].el_id_new,
						&flag_in, &d0, &d1, &d2, &d3, &d4);

					//if not in the element bring it back in the element
					if (flag_in == 0)
					{
						backintria(&prtcls[j]);
						check_el2(prtcls[j].xyz_old[0], prtcls[j].xyz_old[1], prtcls[j].xyz_old[2], prtcls[j].el_id_new,
							&flag_in, &d0, &d1, &d2, &d3, &d4);
					}

					if (prtcls[j].tria_crossed[0] >= 0 && prtcls[j].tria_crossed[1] >= 0
						&& prtcls[j].tria_crossed[2] >= 0 && prtcls[j].node_alone >= 0 && flag_in >= 0)
					{
						//find face that is negative (boundary tria) and return normal
						normal2tria(abs(el_t), j, prtcls[j].v_new[0], prtcls[j].v_new[1], prtcls[j].v_new[2],
							&nor_face_x, &nor_face_y, &nor_face_z);

						//projected velocity vector
						proj_vec(&nor_face_x, &nor_face_y, &nor_face_z, &proj_v_x, &proj_v_y, &proj_v_z,
							prtcls[j].v_new[0], prtcls[j].v_new[1], prtcls[j].v_new[2]);

						//calculate new position
						prtcls[j].xyz_new[0] = prtcls[j].xyz_old[0] + proj_v_x*(t)+nor_face_x*0.001;
						prtcls[j].xyz_new[1] = prtcls[j].xyz_old[1] + proj_v_y*(t)+nor_face_y*0.001;
						prtcls[j].xyz_new[2] = prtcls[j].xyz_old[2] + proj_v_z*(t)+nor_face_z*0.001;

						int el_t_old = el_t;

						//find new element
						el_t = get_block_el(&prtcls[j], &crossed_n1, &crossed_n2, &crossed_n3, &alone_n);

						prtcls[j].tria_crossed[0] = crossed_n1;
						prtcls[j].tria_crossed[1] = crossed_n2;
						prtcls[j].tria_crossed[2] = crossed_n3;
						prtcls[j].node_alone = alone_n;

						//if the particle is still outside the boundaries correct the position again
						if (el_t < 0 && prtcls[j].tria_crossed[0] >= 0 && prtcls[j].tria_crossed[1] >= 0
							&& prtcls[j].tria_crossed[2] >= 0 && prtcls[j].node_alone >= 0)
						{
							//find face that is negative (boundary tria) and return normal
							normal2tria(abs(el_t), j, prtcls[j].v_new[0], prtcls[j].v_new[1], prtcls[j].v_new[2],
								&nor_face_x, &nor_face_y, &nor_face_z);

							//projected velocity vector
							proj_vec(&nor_face_x, &nor_face_y, &nor_face_z, &proj_v_x, &proj_v_y, &proj_v_z,
								prtcls[j].v_new[0], prtcls[j].v_new[1], prtcls[j].v_new[2]);

							prtcls[j].xyz_new[0] = prtcls[j].xyz_old[0] + proj_v_x*(t)+nor_face_x*0.001;
							prtcls[j].xyz_new[1] = prtcls[j].xyz_old[1] + proj_v_y*(t)+nor_face_y*0.001;
							prtcls[j].xyz_new[2] = prtcls[j].xyz_old[2] + proj_v_z*(t)+nor_face_z*0.001;

							el_t = get_block_el(&prtcls[j], &crossed_n1, &crossed_n2, &crossed_n3, &alone_n);

							prtcls[j].tria_crossed[0] = crossed_n1;
							prtcls[j].tria_crossed[1] = crossed_n2;
							prtcls[j].tria_crossed[2] = crossed_n3;
							prtcls[j].node_alone = alone_n;
						}
					}

					if (el_t < 0)
					{
						prtcls[j].xyz_new[0] = temp_x;
						prtcls[j].xyz_new[1] = temp_y;
						prtcls[j].xyz_new[2] = temp_z;
						prtcls[j].track = 0;
					}
				}

				if (el_t > 0)
				{
					//copy new_el into old_el
					prtcls[j].el_id_old = prtcls[j].el_id_new;

					//assign new landing element
					prtcls[j].el_id_new = el_t;

					//copy new barycentric coordinates into old
					for (i = 0; i < 4; i++)
						prtcls[j].bcoord_old[i] = prtcls[j].bcoord_new[i];

					//calculate new barycentric coordinates-> use baricentric
					baricentric(&prtcls[j]);

					//calculate new velocity-> use avg_v
					avg_v(&prtcls[j]);

					charge_els[el_t].id = 1;

				}
			}
		}
	}

	//return 1;
}

EXTERN void write_full_charge(int type)
{
	FILE *Charge_Weld_Full;
	int i = 0;
	if (type == 1)
		Charge_Weld_Full = fopen("Charge_Weld_Full.vtk", "w");
	else
		Charge_Weld_Full = fopen("Seam_Welds_Full.vtk", "w");

	fprintf(Charge_Weld_Full, "# vtk DataFile Version 2.0");
	fprintf(Charge_Weld_Full, "\nExtrusion Processtime %f", nl_Time);
	fprintf(Charge_Weld_Full, "\nASCII");
	fprintf(Charge_Weld_Full, "\n\nDATASET UNSTRUCTURED_GRID");
	fprintf(Charge_Weld_Full, "\nPOINTS %i float", me_nNodes);
	for (i = 0; i<me_nNodes; i++)
	{
		fprintf(Charge_Weld_Full, "\n %f %f %f", me_pNdata[i].Coord[0], me_pNdata[i].Coord[1], me_pNdata[i].Coord[2]);
	}
	fprintf(Charge_Weld_Full, "\n\nCELLS %i %d", me_nTetras, 5 * me_nTetras);
	for (i = 0; i<me_nTetras; i++)
	{
		fprintf(Charge_Weld_Full, "\n4 %i %i %i %i", me_pETdata[i].Nodes[0], me_pETdata[i].Nodes[1], me_pETdata[i].Nodes[2], me_pETdata[i].Nodes[3]);
	}
	fprintf(Charge_Weld_Full, "\n\nCELL_TYPES %i", me_nTetras);
	for (i = 0; i<me_nTetras; i++)
	{
		fprintf(Charge_Weld_Full, "\n10");
	}
	fprintf(Charge_Weld_Full, "\nCELL_DATA %i", me_nTetras);
	fprintf(Charge_Weld_Full, "\nSCALARS Charge float 1");
	fprintf(Charge_Weld_Full, "\nLOOKUP_TABLE default");
	for (i = 0; i<me_nTetras; i++)
	{

		fprintf(Charge_Weld_Full, "\n%i", charge_els[i].id);
	}
	fclose(Charge_Weld_Full);
}

EXTERN void check_trias_area()
{
	int i = 0;
	double area = 0;
	A_mesh = 0;
	for (i = 0; i<rm_nTTrias; i++)
	{
		if (rm_pTTdata[i].Area_cons == 1)
		{
			tria_area(i, &area);
			A_mesh += area;
			if (A_mesh > A_true) A_mesh = A_true;
		}
	}
}

void connectivity(int plane_el[], int el_act_regions[], int el_regions[], int num_el, int *num_holes)
{
	int i = 0, j = 0, connectivity_index = 1, neighbours[4], index[4], added = 0, start = 0, counter = 0, search_index = 0;
	int search = 1;

	//start regions search
	while (search)
	{
		added = 0;
		counter = 0;
		if (search_index == 0)
		{
			el_act_regions[search_index] = 1;
			el_regions[search_index] = connectivity_index;
		}
		else
		{
			el_act_regions[search_index] = 1;
		}
		//get neighbours
		for (j = 0; j<4; j++)
		{
			neighbours[j] = me_pETdata[plane_el[search_index]].Nbr[j];
		}
		//check if the neighbors are part of plane_el
		for (j = 0; j<4; j++)
		{
			if (neighbours[j] > 0)
			{
				for (i = 0; i<num_el; i++)
				{
					if (neighbours[j] == plane_el[i]) break;
				}
				if (i == num_el) index[j] = 10000000; else index[j] = i;
			}
			else
				index[j] = 10000000;
		}
		//check if neighbors have been already added to the el_regions
		// if not add them 
		for (j = 0; j<4; j++)
		{
			if (index[j] < 10000000)
			{
				if (el_act_regions[index[j]] == 0)
				{
					added++;
					el_regions[index[j]] = connectivity_index;
				}
			}
		}
		for (i = 0; i<num_el; i++)
		{
			if (el_regions[i] == connectivity_index && el_act_regions[i] == 0)
				added++;
		}
		//change connectivity index if no neighbors have been added
		if (added == 0)
		{
			connectivity_index++;
			for (i = 0; i<num_el; i++)
			{
				if (el_regions[i] == 0)
				{
					el_regions[i] = connectivity_index;
					search_index = i;
					break;
				}
			}
		}
		else
		{
			for (i = 0; i<num_el; i++)
			{
				if (el_regions[i] == connectivity_index && el_act_regions[i] == 0)
				{
					search_index = i;
					break;
				}
			}
		}

		//break criterium, count for 0s in el_regions
		for (i = 0; i<num_el; i++)
			if (el_regions[i] == 0)
				counter++;
		if (counter == 0) search = 0;
	}
	*num_holes = connectivity_index;

	if (nl_WriteSWConnectivityToCSV) {
		FILE *debug_connectivity;

		printf("\n po: The connectivity search of the portholes has terminated");

		debug_connectivity = fopen("Connectivity_debug.vtk", "w");
		fprintf(debug_connectivity, "# vtk DataFile Version 2.0");
		fprintf(debug_connectivity, "\nExtrusion Processtime %f", nl_Time);
		fprintf(debug_connectivity, "\nASCII");
		fprintf(debug_connectivity, "\n\nDATASET UNSTRUCTURED_GRID");
		fprintf(debug_connectivity, "\nPOINTS %i float", me_nNodes);

		for (i = 0; i < me_nNodes; i++)
		{
			fprintf(debug_connectivity, "\n %f %f %f", me_pNdata[i].Coord[0], me_pNdata[i].Coord[1], me_pNdata[i].Coord[2]);
		}
		fprintf(debug_connectivity, "\nCELLS %i %d", me_nTetras, 5 * me_nTetras);
		for (i = 0; i < me_nTetras; i++)
		{
			fprintf(debug_connectivity, "\n4 %i %i %i %i", me_pETdata[i].Nodes[0], me_pETdata[i].Nodes[1], me_pETdata[i].Nodes[2], me_pETdata[i].Nodes[3]);
		}
		fprintf(debug_connectivity, "\n\nCELL_TYPES %i", me_nTetras);
		for (i = 0; i < me_nTetras; i++)
		{
			fprintf(debug_connectivity, "\n10");
		}
		fprintf(debug_connectivity, "\n\nCELL_DATA %i", me_nTetras);
		fprintf(debug_connectivity, "\nSCALARS Connectivity float 1");
		fprintf(debug_connectivity, "\nLOOKUP_TABLE default");
		for (i = 0; i < me_nTetras; i++)
		{
			for (j = 0; j < num_el; j++)
			{
				if (i == plane_el[j])
				{
					fprintf(debug_connectivity, "\n%i", el_regions[j]);
					break;
				}
			}
			if (j == num_el)
				fprintf(debug_connectivity, "\n%i", 0);
		}
		fclose(debug_connectivity);
	}
}

void write_exit_sw_id(int n_iter)
{
	FILE *pFile;

	int i,j;

	int *n_particles_node;
	double *SW_nodes_qual_Q, *SW_nodes_qual_K, *SW_nodes_qual_J;
	double *SW_nodes_qual_Q_max, *SW_nodes_qual_K_max, *SW_nodes_qual_J_max;

	
	
	n_particles_node = (int*)malloc(sizeof(int)*rm_nTNodes);

	SW_nodes_qual_Q = (double*)malloc(sizeof(double)*rm_nTNodes);
	SW_nodes_qual_K = (double*)malloc(sizeof(double)*rm_nTNodes);
	SW_nodes_qual_J = (double*)malloc(sizeof(double)*rm_nTNodes);
	SW_nodes_qual_Q_max = (double*)malloc(sizeof(double)*rm_nTNodes);
	SW_nodes_qual_K_max = (double*)malloc(sizeof(double)*rm_nTNodes);
	SW_nodes_qual_J_max = (double*)malloc(sizeof(double)*rm_nTNodes);

	for (i = 0; i < rm_nTNodes; i++)
	{
		n_particles_node[i] = 0;
		SW_nodes_qual_Q[i] = 0;
		SW_nodes_qual_K[i] = 0;
		SW_nodes_qual_J[i] = 0;
		SW_nodes_qual_Q_max[i] = 0;
		SW_nodes_qual_K_max[i] = 0;
		SW_nodes_qual_J_max[i] = 0;
	}


	for (j = 0; j < n_prtcls; j++)
	{
		if (prtcls[j].sw_cons > 0)
		{
			int n1 = -1, n2 = -1, n3 = -1, elem = -1, nd = -1, node = -10;
			double l1 = -1, l2 = -1, l3 = -1;

			elem = prtcls[j].el_sw;

			n1 = rm_pNodeNew2Old[rm_pTTdata[elem].TNodes[0]];
			n2 = rm_pNodeNew2Old[rm_pTTdata[elem].TNodes[1]];
			n3 = rm_pNodeNew2Old[rm_pTTdata[elem].TNodes[2]];

			l1 = DIST3(me_pNdata[n1].Coord[0], me_pNdata[n1].Coord[1], me_pNdata[n1].Coord[2],
				prtcls[j].xyz_new[0], prtcls[j].xyz_new[1], prtcls[j].xyz_new[2]);
			l2 = DIST3(me_pNdata[n2].Coord[0], me_pNdata[n2].Coord[1], me_pNdata[n2].Coord[2],
				prtcls[j].xyz_new[0], prtcls[j].xyz_new[1], prtcls[j].xyz_new[2]);
			l3 = DIST3(me_pNdata[n3].Coord[0], me_pNdata[n3].Coord[1], me_pNdata[n3].Coord[2],
				prtcls[j].xyz_new[0], prtcls[j].xyz_new[1], prtcls[j].xyz_new[2]);

			if (l1 < l2 && l1 < l3)
				node = rm_pTTdata[elem].TNodes[0];
			else if (l2 < l1 && l2 < l3)
				node = rm_pTTdata[elem].TNodes[1];
			else if (l3 < l1 && l3 < l2)
				node = rm_pTTdata[elem].TNodes[2];

			n_particles_node[node] = n_particles_node[node]+1;
			SW_nodes_qual_Q[node] += prtcls[j].qual_sw_Q;
			SW_nodes_qual_K[node] += prtcls[j].qual_sw_K;
			SW_nodes_qual_J[node] += prtcls[j].qual_sw_J;

			if (SW_nodes_qual_Q_max[node] < prtcls[j].qual_sw_Q)
			{
				SW_nodes_qual_Q_max[node] = prtcls[j].qual_sw_Q;
			}
			if (SW_nodes_qual_K_max[node] < prtcls[j].qual_sw_K)
			{
				SW_nodes_qual_K_max[node] = prtcls[j].qual_sw_K;
			}
			if (SW_nodes_qual_J_max[node] < prtcls[j].qual_sw_J)
			{
				SW_nodes_qual_J_max[node] = prtcls[j].qual_sw_J;
			}
		}
	}

	for (i = 0; i < rm_nTNodes; i++)
	{
		if (n_particles_node[i] > 0)
		{
			SW_nodes_qual_Q[i] /= n_particles_node[i];
			SW_nodes_qual_K[i] /= n_particles_node[i];
			SW_nodes_qual_J[i] /= n_particles_node[i];
		}
	}

	char Name[255];
	sprintf(Name, "SW_Exit_qual_%i.vtk", n_iter);

	pFile = fopen(Name, "w");

	fprintf(pFile, "%s\n", "# vtk DataFile Version 1.0");
	fprintf(pFile, "%s\n", "Unstructured Grid Example");
	fprintf(pFile, "%s\n\n", "ASCII");
	fprintf(pFile, "%s\n", "DATASET UNSTRUCTURED_GRID");
	fprintf(pFile, "%s  %d  float\n", "POINTS", rm_nTNodes);

	
	for (i = 0; i < rm_nTNodes; i++)
	{
		int nd = rm_pNodeNew2Old[i];
		fprintf(pFile, "%10.3f  %10.3f  %10.3f\n",
			me_pNdata[nd].Coord[0], me_pNdata[nd].Coord[1], me_pNdata[nd].Coord[2]);
	}
	fprintf(pFile, "\n");

	fprintf(pFile, "%s %d  %d\n", "CELLS", rm_nTTrias, rm_nTTrias * 4);

	int e;
	for (e = 0; e < rm_nTTrias; e++)
		fprintf(pFile, "3 %8d  %8d  %8d\n",
			rm_pTTdata[e].TNodes[0], rm_pTTdata[e].TNodes[1],
			rm_pTTdata[e].TNodes[2]);
	fprintf(pFile, "\n");

	fprintf(pFile, "%s %d\n", "CELL_TYPES", rm_nTTrias);
	for (e = 0; e < rm_nTTrias; e++)
		fprintf(pFile, "5 \n");
	fprintf(pFile, "\n\n");

	fprintf(pFile, "POINT_DATA %d\n", rm_nTNodes);
	fprintf(pFile, "SCALARS Q-Value-avg float\n");
	fprintf(pFile, "LOOKUP_TABLE default\n");
	for (i = 0; i < rm_nTNodes; i++)
	{
		fprintf(pFile, "%f \n", SW_nodes_qual_Q[i]);
	}

	fprintf(pFile, "SCALARS K-Value-avg float\n");
	fprintf(pFile, "LOOKUP_TABLE default\n");
	for (i = 0; i < rm_nTNodes; i++)
	{
		fprintf(pFile, "%f \n", SW_nodes_qual_K[i]);
	}

	fprintf(pFile, "SCALARS J-Value-avg float\n");
	fprintf(pFile, "LOOKUP_TABLE default\n");
	for (i = 0; i < rm_nTNodes; i++)
	{
		fprintf(pFile, "%f \n", SW_nodes_qual_J[i]);
	}

	fprintf(pFile, "SCALARS Q-Value-max float\n");
	fprintf(pFile, "LOOKUP_TABLE default\n");
	for (i = 0; i < rm_nTNodes; i++)
	{
		fprintf(pFile, "%f \n", SW_nodes_qual_Q_max[i]);
	}

	fprintf(pFile, "SCALARS K-Value-max float\n");
	fprintf(pFile, "LOOKUP_TABLE default\n");
	for (i = 0; i < rm_nTNodes; i++)
	{
		fprintf(pFile, "%f \n", SW_nodes_qual_K_max[i]);
	}

	fprintf(pFile, "SCALARS J-Value-max float\n");
	fprintf(pFile, "LOOKUP_TABLE default\n");
	for (i = 0; i < rm_nTNodes; i++)
	{
		fprintf(pFile, "%f \n", SW_nodes_qual_J_max[i]);
	}

	fclose(pFile);

}

void sort_exit_sw_id()
{
	int i, j;
	for (i = 0; i < rm_nTTrias; i++)
	{
		if (id_TTrias_parsing[i] != 0)
		{
			int counter = 0;
			if (id_TTrias_parsing[i] < 1000)
				for (j = 0; j < 3; j++)
					if (rm_pTTdata[i].NeighborTria[j] < rm_nTTrias)
						if (id_TTrias_parsing[rm_pTTdata[i].NeighborTria[j]] == id_TTrias_parsing[i])
							counter++;
			if (counter > 2)
				id_TTrias_parsing[i] = 0;
			else
			{
				if (id_TTrias_parsing[i] > 1000)
				{
					int n1 = rm_pTTdata[i].Nodes[0];
					int n2 = rm_pTTdata[i].Nodes[1];
					int n3 = rm_pTTdata[i].Nodes[2];
					int n1t = rm_pTTdata[i].TNodes[0];
					int n2t = rm_pTTdata[i].TNodes[1];
					int n3t = rm_pTTdata[i].TNodes[2];

					double n1x = me_pNdata[n1].Coord[0]; double n1y = me_pNdata[n1].Coord[1];

					double n2x = me_pNdata[n2].Coord[0]; double n2y = me_pNdata[n2].Coord[1];

					double n3x = me_pNdata[n3].Coord[0]; double n3y = me_pNdata[n3].Coord[1];

					double r1_1 = DIST2(n1x, n1y, n2x, n2y) / 3;

					double r1_2 = DIST2(n1x, n1y, n3x, n3y) / 3;

					double r2_1 = DIST2(n2x, n2y, n1x, n1y) / 3;

					double r2_2 = DIST2(n2x, n2y, n3x, n3y) / 3;

					double r3_1 = DIST2(n3x, n3y, n1x, n1y) / 3;

					double r3_2 = DIST2(n3x, n3y, n2x, n2y) / 3;

					double r1 = 0, r2 = 0, r3 = 0;

					if (r1_1 < r1_2)
						r1 = r1_1;
					else
						r1 = r1_2;

					if (r2_1 < r2_2)
						r2 = r2_1;
					else
						r2 = r2_2;

					if (r3_1 < r3_2)
						r3 = r3_1;
					else
						r3 = r3_2;

					for (j = 0; j < n_prtcls; j++)
					{
						if (prtcls[j].track != 0)
						{
							if (prtcls[j].el_sw == i)
							{
								double px = prtcls[j].xyz_new[0];
								double py = prtcls[j].xyz_new[1];
								double test = DIST2(n1x, n1y, px, py);
								if (DIST2(n1x, n1y, px, py) > r1)
									if (id_Tnodes_parsing[n1t] == 0)
										id_Tnodes_parsing[n1t] = prtcls[j].id_sw;
									else
										id_Tnodes_parsing[n1t] = id_Tnodes_parsing[n1t] * 1000 + prtcls[j].id_sw;
								else if (DIST2(n2x, n2y, px, py) > r2)
									if (id_Tnodes_parsing[n2t] == 0)
										id_Tnodes_parsing[n2t] = prtcls[j].id_sw;
									else
										id_Tnodes_parsing[n2t] = id_Tnodes_parsing[n2t] * 1000 + prtcls[j].id_sw;
								else if (DIST2(n3x, n3y, px, py) > r3)
									if (id_Tnodes_parsing[n3t] == 0)
										id_Tnodes_parsing[n3t] = prtcls[j].id_sw;
									else
										id_Tnodes_parsing[n3t] = id_Tnodes_parsing[n3t] * 1000 + prtcls[j].id_sw;
							}
						}
					}
				}
				else
				{
					counter = 0;
					int el = -1;
					if (rm_pTTdata[i].NeighborTria[0] < rm_nTTrias)
						if (id_TTrias_parsing[rm_pTTdata[i].NeighborTria[0]] != id_TTrias_parsing[i] && id_TTrias_parsing[rm_pTTdata[i].NeighborTria[0]] < 1000
							&& id_TTrias_parsing[rm_pTTdata[i].NeighborTria[0]] > 0)
						{
							counter++;
							el = 0;
						}
					if (rm_pTTdata[i].NeighborTria[1] < rm_nTTrias)
						if (id_TTrias_parsing[rm_pTTdata[i].NeighborTria[1]] != id_TTrias_parsing[i] && id_TTrias_parsing[rm_pTTdata[i].NeighborTria[1]] < 1000
							&& id_TTrias_parsing[rm_pTTdata[i].NeighborTria[1]] > 0)
						{
							counter++;
							el = 1;
						}
					if (rm_pTTdata[i].NeighborTria[2] < rm_nTTrias)
						if (id_TTrias_parsing[rm_pTTdata[i].NeighborTria[2]] != id_TTrias_parsing[i] && id_TTrias_parsing[rm_pTTdata[i].NeighborTria[2]] < 1000
							&& id_TTrias_parsing[rm_pTTdata[i].NeighborTria[2]] > 0)
						{
							counter++;
							el = 2;
						}
					if (el > -1)
					{
						int m, n;
						for (m = 0; m < 3; m++)
						{
							int node = rm_pTTdata[i].Nodes[m];
							/*						node = rm_pNodeNew2Old[node];*/
							for (n = 0; n < 3; n++)
							{
								int node2 = rm_pTTdata[rm_pTTdata[i].NeighborTria[el]].Nodes[n];
								//node2 = rm_pNodeNew2Old[node2];
								if (node == node2)
								{
									int nd = rm_pTTdata[i].TNodes[m];
									id_Tnodes_parsing[nd] = id_TTrias_parsing[i] * 1000 + id_TTrias_parsing[rm_pTTdata[i].NeighborTria[el]];
								}
							}
						}
					}
				}
			}
		}
	}
	FILE *SW;
	SW = fopen("SW_Points_DEBUG.csv", "w");
	fprintf(SW, "x,y,z");
	for (j = 0; j < rm_nTNodes; j++)
		if (id_Tnodes_parsing[j] > 1000)
		{
			int nd = rm_pNodeNew2Old[j];
			fprintf(SW, "\n%10.3f,  %10.3f,  %10.3f\n", me_pNdata[nd].Coord[0], me_pNdata[nd].Coord[1], me_pNdata[nd].Coord[2]);
		}
	fclose(SW);
}


EXTERN int am_po_bearingopt()
{
	int nFoundNodes = 0;
	int regionID = -1;
	double tol = 0.20;
	int debug = 1;        // activate extra output for debugging purposes
	int bandwidth = 0;    // activate bandwidth mode
	int terminate = 1;    // terminate optimization

	nl_counterAndrea++;
	bearingInfos *BearingOptInfo;


	// 1. Find pBearingRegionFlag index
	for (int i = 0; i < me_nTriaReg; i++) {
		if (nl_pPdata[0].pBearingRegionFlag[i] == 1) regionID = i;
	}

	// 2. Find z_max
	double z_max = me_pNdata[me_pTriaReg[0].pTriaNodes[0][0]].OldCoord[2];
	if (regionID > -1)
	{
		if (nl_PunchVelocity > 0) {
			for (int m = 0; m < me_nTriaReg; m++) {
				if (nl_pPdata[0].pBearingRegionFlag[m] == 1) {
					for (int i = 0; i < me_pTriaReg[m].nTrias; i++) {
						for (int k = 0; k < 3; k++) {
							if (me_pNdata[me_pTriaReg[m].pTriaNodes[i][k]].OldCoord[2] > z_max) {
								z_max = me_pNdata[me_pTriaReg[m].pTriaNodes[i][k]].OldCoord[2];
							}
						}
					}
				}
			}
		}
		else if (nl_PunchVelocity < 0) {
			for (int m = 0; m < me_nTriaReg; m++) {
				if (nl_pPdata[0].pBearingRegionFlag[m] == 1) {
					for (int i = 0; i < me_pTriaReg[m].nTrias; i++) {
						for (int k = 0; k < 3; k++) {
							if (me_pNdata[me_pTriaReg[m].pTriaNodes[i][k]].OldCoord[2] < z_max) {
								z_max = me_pNdata[me_pTriaReg[m].pTriaNodes[i][k]].OldCoord[2];
							}
						}
					}
				}
			}

		}
	}
	else // in case no pBearingRegionFlag index present
	{
		if (nl_PunchVelocity > 0) {
			for (int m = 0; m < me_nTriaReg; m++) {
				for (int i = 0; i < me_pTriaReg[m].nTrias; i++) {
					for (int k = 0; k < 3; k++) {
						if (me_pNdata[me_pTriaReg[m].pTriaNodes[i][k]].OldCoord[2] > z_max) {
							z_max = me_pNdata[me_pTriaReg[m].pTriaNodes[i][k]].OldCoord[2];
							nl_pPdata[0].pBearingRegionFlag[m] = 1;
							regionID++;
						}
					}
				}
			}
		}
		else if (nl_PunchVelocity < 0) {
			for (int m = 0; m < me_nTriaReg; m++) {
				for (int i = 0; i < me_pTriaReg[m].nTrias; i++) {
					for (int k = 0; k < 3; k++) {
						if (me_pNdata[me_pTriaReg[m].pTriaNodes[i][k]].OldCoord[2] < z_max) {
							z_max = me_pNdata[me_pTriaReg[m].pTriaNodes[i][k]].OldCoord[2];
							nl_pPdata[0].pBearingRegionFlag[m] = 1;
							regionID++;
						}
					}
				}
			}

		}
	}

	// 3. Find z_min
	double z_min = z_max;
	if (regionID > -1)
	{
		if (nl_PunchVelocity > 0) {
			for (int m = 0; m < me_nTriaReg; m++) {
				if (nl_pPdata[0].pBearingRegionFlag[m] == 1) {
					for (int i = 0; i < me_pTriaReg[m].nTrias; i++) {
						for (int k = 0; k < 3; k++) {
							if (me_pNdata[me_pTriaReg[m].pTriaNodes[i][k]].OldCoord[2] < z_min) {
								z_min = me_pNdata[me_pTriaReg[m].pTriaNodes[i][k]].OldCoord[2];
							}
						}
					}
				}
			}
		}
		else if (nl_PunchVelocity < 0) {
			for (int m = 0; m < me_nTriaReg; m++) {
				if (nl_pPdata[0].pBearingRegionFlag[m] == 1) {
					for (int i = 0; i < me_pTriaReg[m].nTrias; i++) {
						for (int k = 0; k < 3; k++) {
							if (me_pNdata[me_pTriaReg[m].pTriaNodes[i][k]].OldCoord[2] > z_min) {
								z_min = me_pNdata[me_pTriaReg[m].pTriaNodes[i][k]].OldCoord[2];
							}
						}
					}
				}
			}

		}
	}

	// 4. Search for bearing line (Nodes that are delimiting the bearing)
	if (regionID > -1)
	{
		for (int m = regionID; m <= regionID; m++) {
			if (nl_pPdata[0].pBearingRegionFlag[m] == 1) {
				for (int i = 0; i < me_pTriaReg[m].nTrias; i++) {
					for (int k = 0; k < 3; k++) {
						int count_shared_node = 0;
						int consider = 1;
						for (int j = regionID; j <= regionID; j++) {
							for (int l = 0; l < me_pTriaReg[j].nTrias; l++) {
								for (int n = 0; n < 3; n++) {
									
									if (me_pTriaReg[m].pTriaNodes[i][k] == me_pTriaReg[j].pTriaNodes[l][n]) {
										count_shared_node = count_shared_node + 1;
										if (me_pTriaReg[j].BearingLine[l][n] != -1) consider = 0;
									}
								}
							}
						}
						

						if (count_shared_node < 6 && me_pNdata[me_pTriaReg[m].pTriaNodes[i][k]].OldCoord[2] > z_min && consider == 1) {
							me_pTriaReg[m].BearingLine[i][k] = me_pTriaReg[m].pTriaNodes[i][k];
							nFoundNodes++;
						}

					}
				}
			}

		}
	}

	
	if (debug == 1) {
		FILE *found_points1;
		char vtkFileName[255];

		sprintf(vtkFileName, "found_points1_%i.csv", nl_counterAndrea);

		found_points1 = fopen(vtkFileName, "w");

		fprintf(found_points1, "x,y,z\n");

		for (int i = 0; i < me_pTriaReg[regionID].nTrias; i++) {
			for (int k = 0; k < 3; k++) {
				if (me_pTriaReg[regionID].BearingLine[i][k] != -1) {
					double x = me_pNdata[me_pTriaReg[regionID].BearingLine[i][k]].OldCoord[0];
					double y = me_pNdata[me_pTriaReg[regionID].BearingLine[i][k]].OldCoord[1];
					double z = me_pNdata[me_pTriaReg[regionID].BearingLine[i][k]].OldCoord[2];
					fprintf(found_points1, "%f,%f,%f\n", x, y, z);
				}
			}
		}
		fclose(found_points1);
	}

	//extra condition for bearing node search
	if (regionID > -1)
	{
		for (int i = 0; i < me_nTriaReg; i++) {
			if (nl_pPdata[0].pBearingRegionFlag[i] == 1) {
				for (int j = 0; j < me_pTriaReg[i].nTrias; j++) {
					int consider = 0;
					int nodo[3];
					nodo[0] = -1;
					nodo[1] = -1;
					nodo[2] = -1;
					for (int k = 0; k < 3; k++) {
						for (int m = 0; m < me_nTriaReg; m++) {
							for (int n = 0; n < me_pTriaReg[m].nTrias; n++) {
								for (int l = 0; l < 3; l++) {
									if (me_pTriaReg[m].BearingLine[n][l] != -1 &&
										me_pTriaReg[m].pTriaNodes[n][l] == me_pTriaReg[i].pTriaNodes[j][k]) {
										nodo[consider] = k;
										consider++;
									}
								}
							}
						}
					}
					
					if (consider == 3) {
						nFoundNodes--;
						double z_min = me_pNdata[me_pTriaReg[i].pTriaNodes[j][0]].OldCoord[2];
						int rem = me_pTriaReg[i].pTriaNodes[j][0];
						for (int a = 0; a < 3; a++) {
							double z = me_pNdata[me_pTriaReg[i].pTriaNodes[j][a]].OldCoord[2];
							if (z < z_min) {
								z_min = z;
								rem = me_pTriaReg[i].pTriaNodes[j][a];
							}
						}
						for (int c = 0; c < me_pTriaReg[i].nTrias; c++) {
							for (int a = 0; a < 3; a++) {
								if (me_pTriaReg[i].pTriaNodes[c][a] == rem) {
									me_pTriaReg[i].BearingLine[c][a] = -1;
								}
							}
						}
					}
					
					if (consider == 2) {
						me_pTriaReg[i].BearingElement[j] = 1;
						double x1 = me_pNdata[me_pTriaReg[i].pTriaNodes[j][nodo[0]]].OldCoord[0];
						double y1 = me_pNdata[me_pTriaReg[i].pTriaNodes[j][nodo[0]]].OldCoord[1];
						double z1 = me_pNdata[me_pTriaReg[i].pTriaNodes[j][nodo[0]]].OldCoord[2];
						double x2 = me_pNdata[me_pTriaReg[i].pTriaNodes[j][nodo[1]]].OldCoord[0];
						double y2 = me_pNdata[me_pTriaReg[i].pTriaNodes[j][nodo[1]]].OldCoord[1];
						double z2 = me_pNdata[me_pTriaReg[i].pTriaNodes[j][nodo[1]]].OldCoord[2];

						double dist12 = sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));

						for (int b = 0; b < me_pTriaReg[i].nTrias; b++) {
							for (int c = 0; c < 3; c++) {
								double xB = me_pNdata[me_pTriaReg[i].pTriaNodes[b][c]].OldCoord[0];
								double yB = me_pNdata[me_pTriaReg[i].pTriaNodes[b][c]].OldCoord[1];
								double zB = me_pNdata[me_pTriaReg[i].pTriaNodes[b][c]].OldCoord[2];
								if (me_pTriaReg[i].BearingLine[b][c] != -1 && zB < z2 && zB < z1) {
									double dist1 = sqrt(pow(x1 - xB, 2) + pow(y1 - yB, 2));
									double dist2 = sqrt(pow(x2 - xB, 2) + pow(y2 - yB, 2));

									if (dist1 <= dist12 && dist2 <= dist12) {
										me_pTriaReg[i].BearingLine[b][c] = -1;
										nFoundNodes--;
									}
								}
							}
						}

					}
				}
				for (int j = 0; j < me_pTriaReg[i].nTrias; j++) {
					for (int k = 0; k < 3; k++) {
						if (me_pTriaReg[i].BearingLine[j][k] != -1) {
							double x1 = me_pNdata[me_pTriaReg[i].pTriaNodes[j][k]].OldCoord[0];
							double y1 = me_pNdata[me_pTriaReg[i].pTriaNodes[j][k]].OldCoord[1];
							double z1 = me_pNdata[me_pTriaReg[i].pTriaNodes[j][k]].OldCoord[2];
							for (int a = 0; a < me_pTriaReg[i].nTrias; a++) {
								for (int b = 0; b < 3; b++) {
									if (me_pTriaReg[i].BearingLine[a][b] != -1) {
										double x2 = me_pNdata[me_pTriaReg[i].pTriaNodes[a][b]].OldCoord[0];
										double y2 = me_pNdata[me_pTriaReg[i].pTriaNodes[a][b]].OldCoord[1];
										double z2 = me_pNdata[me_pTriaReg[i].pTriaNodes[a][b]].OldCoord[2];
										if (x1 == x2 && y1 == y2 && z1 < z2) {
											me_pTriaReg[i].BearingLine[j][k] = -1;
											nFoundNodes--;
										}
									}
								}
							}

						}
					}
				}
			}
		}
	}


	if (debug == 1) {
		FILE *found_points;
		char vtkFileName[255];

		sprintf(vtkFileName, "found_points_%i.csv", nl_counterAndrea);

		found_points = fopen(vtkFileName, "w");

		fprintf(found_points, "x,y,z\n");

		for (int i = 0; i < me_pTriaReg[regionID].nTrias; i++) {
			for (int k = 0; k < 3; k++) {
				if (me_pTriaReg[regionID].BearingLine[i][k] != -1) {
					double x = me_pNdata[me_pTriaReg[regionID].BearingLine[i][k]].OldCoord[0];
					double y = me_pNdata[me_pTriaReg[regionID].BearingLine[i][k]].OldCoord[1];
					double z = me_pNdata[me_pTriaReg[regionID].BearingLine[i][k]].OldCoord[2];
					fprintf(found_points, "%f,%f,%f\n", x, y, z);
				}
			}
		}
		fclose(found_points);
	}


	// 5. Initilize BearingOpt Struct for information storing
	if ((BearingOptInfo = (bearingInfos *)GE_MALLOC(bearingInfos, [1], nFoundNodes)) == NULL)
		return nl_Errors(nl_AllocError, nFoundNodes * sizeof(bearingInfos), NULL, NULL, NULL);

	
	int l = 0;
	if (regionID > -1)
	{
		for (int m = 0; m < me_nTriaReg; m++) {
			if (nl_pPdata[0].pBearingRegionFlag[m] == 1) {
				for (int i = 0; i < me_pTriaReg[m].nTrias; i++) {
					for (int k = 0; k < 3; k++) {
						if (me_pTriaReg[m].BearingLine[i][k] != -1) {
							BearingOptInfo[l].ID = me_pTriaReg[m].BearingLine[i][k];
							BearingOptInfo[l].BearNodesXYZ[0] = me_pNdata[me_pTriaReg[m].BearingLine[i][k]].OldCoord[0];
							BearingOptInfo[l].BearNodesXYZ[1] = me_pNdata[me_pTriaReg[m].BearingLine[i][k]].OldCoord[1];
							BearingOptInfo[l].BearNodesXYZ[2] = me_pNdata[me_pTriaReg[m].BearingLine[i][k]].OldCoord[2];
							BearingOptInfo[l].ClosestBearingNode = -1;
							l++;
						}
					}
				}
			}
		}
	}


	// 6. Calculation of relative velocity
	double R = nl_ExtrusionRatio;

	for (int i = 0; i < rm_nTTrias; i++) {
		for (int k = 0; k < 3; k++) {
			double v1 = me_pNdata[rm_pTTdata[i].Nodes[k]].Velocity[0];
			double v2 = me_pNdata[rm_pTTdata[i].Nodes[k]].Velocity[1];
			double v3 = me_pNdata[rm_pTTdata[i].Nodes[k]].Velocity[2];
			double v = sqrt(v3*v3);
			rm_pTTdata[i].RelVelocity[k] = v / (R*nl_PunchVelocity);
		}
	}

	if (debug == 1) {
		FILE *relative_velocities;
		char vtkFileName[255];

		sprintf(vtkFileName, "relative_velocities_%i.csv", nl_counterAndrea);

		relative_velocities = fopen(vtkFileName, "w");

		fprintf(relative_velocities, "x,y,z,v_r\n");

		for (int i = 0; i < rm_nTTrias; i++) {
			for (int k = 0; k < 3; k++) {
				double x = me_pNdata[rm_pTTdata[i].Nodes[k]].OldCoord[0];
				double y = me_pNdata[rm_pTTdata[i].Nodes[k]].OldCoord[1];
				double z = me_pNdata[rm_pTTdata[i].Nodes[k]].OldCoord[2];
				double v_r = rm_pTTdata[i].RelVelocity[k];
				fprintf(relative_velocities, "%f,%f,%f,%f\n", x, y, z, v_r);
			}
		}
		fclose(relative_velocities);
	}

	// 7. Referencing of section Nodes to bearing nodes
	for (int i = 0; i < rm_nTTrias; i++) {
		for (int k = 0; k < 3; k++) {
			int checked = 0;
			for (int l = 0; l < rm_nTTrias; l++) {
				for (int n = 0; n < 3; n++) {
					if (rm_pTTdata[i].Nodes[k] == rm_pTTdata[l].Nodes[n] && rm_pTTdata[l].checked[n] == 1) {
						checked = checked + 1;
					}
				}
			}
			if (checked == 0) {
				double dx = me_pNdata[rm_pTTdata[i].Nodes[k]].OldCoord[0] - BearingOptInfo[0].BearNodesXYZ[0];
				double dy = me_pNdata[rm_pTTdata[i].Nodes[k]].OldCoord[1] - BearingOptInfo[0].BearNodesXYZ[1];
				rm_pTTdata[i].BearingDistance[k] = sqrt(dx * dx + dy * dy);
				rm_pTTdata[i].checked[k] = 1;
				for (int j = 0; j < nFoundNodes; j++) {
					double dx = me_pNdata[rm_pTTdata[i].Nodes[k]].OldCoord[0] - BearingOptInfo[j].BearNodesXYZ[0];
					double dy = me_pNdata[rm_pTTdata[i].Nodes[k]].OldCoord[1] - BearingOptInfo[j].BearNodesXYZ[1];
					double dist = sqrt(dx * dx + dy * dy);
					if (dist <= rm_pTTdata[i].BearingDistance[k]) {
						rm_pTTdata[i].BearingDistance[k] = dist;
						rm_pTTdata[i].ClosestBearingNode[k] = BearingOptInfo[j].ID;
					}
				}
			}
		}
	}

	if (debug == 1) {
		FILE *referenced_nodes;
		char vtkFileName[255];

		sprintf(vtkFileName, "referenced_nodes_%i.csv", nl_counterAndrea);

		referenced_nodes = fopen(vtkFileName, "w");

		fprintf(referenced_nodes, "x,y,z,x_b,y_b,z_b\n");
		for (int i = 0; i < rm_nTTrias; i++) {
			for (int k = 0; k < 3; k++) {
				if (rm_pTTdata[i].checked[k] == 1) {
					double x1 = me_pNdata[rm_pTTdata[i].Nodes[k]].OldCoord[0];
					double y1 = me_pNdata[rm_pTTdata[i].Nodes[k]].OldCoord[1];
					double z1 = me_pNdata[rm_pTTdata[i].Nodes[k]].OldCoord[2];
					double x2 = me_pNdata[rm_pTTdata[i].ClosestBearingNode[k]].OldCoord[0];
					double y2 = me_pNdata[rm_pTTdata[i].ClosestBearingNode[k]].OldCoord[1];
					double z2 = me_pNdata[rm_pTTdata[i].ClosestBearingNode[k]].OldCoord[2];
					fprintf(referenced_nodes, "%f,%f,%f,%f,%f,%f\n", x1, y1, z1, x2, y2, z2);
				}
			}
		}
		fclose(referenced_nodes);
	}

	// 8. Calculate equivalent VDR with consideration of distance
	for (int i = 0; i < nFoundNodes; i++) {
		double weight = 0;
		BearingOptInfo[i].BearNodesVDR = 0;
		for (int j = 0; j < rm_nTTrias; j++) {
			for (int k = 0; k < 3; k++) {
				if (rm_pTTdata[j].checked[k] == 1 && rm_pTTdata[j].ClosestBearingNode[k] == BearingOptInfo[i].ID) {
					BearingOptInfo[i].BearNodesVDR = BearingOptInfo[i].BearNodesVDR +
						1 / (1 + rm_pTTdata[j].BearingDistance[k]) * rm_pTTdata[j].RelVelocity[k];
					weight = weight + 1 / (1 + rm_pTTdata[j].BearingDistance[k]);
					BearingOptInfo[i].ClosestBearingNode = 1;
				}
			}
		}
		if (weight != 0) BearingOptInfo[i].BearNodesVDR = BearingOptInfo[i].BearNodesVDR / weight;
	}

	if (debug == 1) {
		FILE *VDR;
		char vtkFileName[255];

		sprintf(vtkFileName, "VDR_%i.csv", nl_counterAndrea);

		VDR = fopen(vtkFileName, "w");

		fprintf(VDR, "x,y,z,VRD\n");
		for (int i = 0; i < nFoundNodes; i++) {
			double x = me_pNdata[BearingOptInfo[i].ID].OldCoord[0];
			double y = me_pNdata[BearingOptInfo[i].ID].OldCoord[1];
			double z = me_pNdata[BearingOptInfo[i].ID].OldCoord[2];
			double VRD = BearingOptInfo[i].BearNodesVDR;
			fprintf(VDR, "%f,%f,%f,%f\n", x, y, z, VRD);
		}
		fclose(VDR);
	}

	// 9. Find VDR_max
	double VDR_max = BearingOptInfo[0].BearNodesVDR;
	for (int i = 0; i < nFoundNodes; i++) {
		if (BearingOptInfo[i].BearNodesVDR > VDR_max) VDR_max = BearingOptInfo[i].BearNodesVDR;
	}

	/*--------------------------------------------------------------------------------------------------------*/
	//                            SHORTENING BEARING LINE   VRD < 1
	/*--------------------------------------------------------------------------------------------------------*/
	int counter = 0;

	if (debug == 1) {
		FILE *Remove_Nodes;
		char vtkFileName[255];

		sprintf(vtkFileName, "Remove Nodes__%i.csv", nl_counterAndrea);

		Remove_Nodes = fopen(vtkFileName, "w");

		fprintf(Remove_Nodes, "x,y,z,\n");
		for (int i = 0; i < nFoundNodes; i++) {
			if (BearingOptInfo[i].BearNodesVDR < 1-tol && BearingOptInfo[i].ClosestBearingNode == 1) {
				double x = me_pNdata[BearingOptInfo[i].ID].OldCoord[0];
				double y = me_pNdata[BearingOptInfo[i].ID].OldCoord[1];
				double z = me_pNdata[BearingOptInfo[i].ID].OldCoord[2];
				fprintf(Remove_Nodes, "%f,%f,%f\n", x, y, z);
			}
		}
		fclose(Remove_Nodes);
	}

	double tol2 = 0.1;

	if (regionID > -1) {
		for (int reg = 0; reg < me_nTriaReg; reg++) {
			if (nl_pPdata[0].pBearingRegionFlag[reg] == 1) {
				int Bearing = -1;
				for (int m = 0; m < co_nTools; m++) {
					if (strcmp(me_pTriaReg[reg].Name, co_pTool[m].ToolName) == 0) Bearing = m;
				}

				for (int i = 0; i < nFoundNodes; i++) {
					double x = me_pNdata[BearingOptInfo[i].ID].OldCoord[0];
					double y = me_pNdata[BearingOptInfo[i].ID].OldCoord[1];
					double z = me_pNdata[BearingOptInfo[i].ID].OldCoord[2];

					//---------------------------------------------------//
					//                     NO BANDWIDTH                  //
					//---------------------------------------------------//
					if (bandwidth == 0) {
						if (BearingOptInfo[i].BearNodesVDR < 1 - tol && BearingOptInfo[i].ClosestBearingNode == 1) {

							// remove all elements containing node
							for (int j = 0; j < co_pTool[Bearing].nFacets; j++) {
								int remove = 0;
								for (int l = 0; l < 3; l++) {
									int A = co_pTool[Bearing].pFacet[j].Points[l];
									double z_A = co_pTool[Bearing].pPCoord[A][2];
									if (z_A > z_min+0.1) remove++;
								}

								if (remove == 3) {              // do not remove last elements
									for (int k = 0; k < 3; k++) {
										int P = co_pTool[Bearing].pFacet[j].Points[k];
										double x_P = co_pTool[Bearing].pPCoord[P][0];
										double y_P = co_pTool[Bearing].pPCoord[P][1];
										double z_P = co_pTool[Bearing].pPCoord[P][2];
										if (fabs(x_P - x) < tol2 && fabs(y_P - y) < tol2 && fabs(z_P - z) < tol2) {
											for (int e = j; e < co_pTool[Bearing].nFacets; e++) {
												if (e != co_pTool[Bearing].nFacets) {
													co_pTool[Bearing].pFacet[e].Points[0] = co_pTool[Bearing].pFacet[e + 1].Points[0];
													co_pTool[Bearing].pFacet[e].Points[1] = co_pTool[Bearing].pFacet[e + 1].Points[1];
													co_pTool[Bearing].pFacet[e].Points[2] = co_pTool[Bearing].pFacet[e + 1].Points[2];
												}
											}
											if ((co_pTool[Bearing].pFacet = (CO_FACET_DATA *)GE_REALLOC(CO_FACET_DATA, [1],
												co_pTool[Bearing].pFacet, co_pTool[Bearing].nFacets - 1)) == NULL)
												return 0; 
											

											co_pTool[Bearing].nFacets -= 1;
											j = 0;
											k = 0;
											terminate = 0;

											//if (!rm_co_ToolInitialisation(12, co_pTool + Bearing)) return 0;
											rm_co_NodeDataInitialisation(me_nBNodes, co_pBNdata);
										}
									}
								}
							}


							// UPDATE me_pTriaReg 
							for (int m = 0; m < me_pTriaReg[reg].nTrias; m++) {
								int remove = 0;
								for (int e = 0; e < 3; e++) {
									double z_A = me_pNdata[me_pTriaReg[reg].pTriaNodes[m][e]].OldCoord[2];
									if (z_A > z_min+0.1) remove++;
								}

								if (remove == 3) {      // do not remove last elements
									for (int n = 0; n < 3; n++) {
										//remove elements containing node
										if (me_pTriaReg[reg].pTriaNodes[m][n] == BearingOptInfo[i].ID) {
											for (int f = m; f < me_pTriaReg[reg].nTrias; f++) {
												if (f != me_pTriaReg[reg].nTrias) {
													me_pTriaReg[reg].pTriaNodes[f][0] = me_pTriaReg[reg].pTriaNodes[f + 1][0];
													me_pTriaReg[reg].pTriaNodes[f][1] = me_pTriaReg[reg].pTriaNodes[f + 1][1];
													me_pTriaReg[reg].pTriaNodes[f][2] = me_pTriaReg[reg].pTriaNodes[f + 1][2];
												}
											}

											if ((me_pTriaReg[reg].pTriaNodes = (int(*)[3])GE_REALLOC(int, [3],
												me_pTriaReg[reg].pTriaNodes, me_pTriaReg[reg].nTrias - 1)) == NULL)
												return 0;

											if ((me_pTriaReg[reg].BearingLine = (int(*)[3])GE_REALLOC(int, [3],
												me_pTriaReg[reg].BearingLine, me_pTriaReg[reg].nTrias - 1)) == NULL)
												return 0;

											if ((me_pTriaReg[reg].BearingElement = (int(*)[1])GE_REALLOC(int, [1],
												me_pTriaReg[reg].BearingElement, me_pTriaReg[reg].nTrias - 1)) == NULL)
												return 0;
												

									
											me_pTriaReg[reg].nTrias -= 1;
											m = 0;
											n = 0;
	
										}
									}
								}
							}
						}
					}
					//---------------------------------------------------//
					//                       BANDWIDTH                   //
					//---------------------------------------------------//
					if (bandwidth == 1) {
						if (BearingOptInfo[i].BearNodesVDR < 0.1) {

							for (int j = 0; j < co_pTool[Bearing].nFacets; j++) {
								int remove = 0;
								for (int l = 0; l < 3; l++) {
									int A = co_pTool[Bearing].pFacet[j].Points[l];
									double z_A = co_pTool[Bearing].pPCoord[A][2];
									if (z_A > z_min+0.1) remove++;
								}
								if (remove == 3) {                                             // avoid removing last elements  //
									for (int k = 0; k < 3; k++) {
										int P = co_pTool[Bearing].pFacet[j].Points[k];
										double x_P = co_pTool[Bearing].pPCoord[P][0];
										double y_P = co_pTool[Bearing].pPCoord[P][1];
										double z_P = co_pTool[Bearing].pPCoord[P][2];
										if (fabs(x_P - x) < tol2 && fabs(y_P - y) < tol2 && fabs(z_P - z) < tol2) {
											int M = co_pTool[Bearing].pFacet[j].Points[0];
											double zM = co_pTool[Bearing].pPCoord[M][2];
											for (int n = 0; n < 3; n++) {
												int N = co_pTool[Bearing].pFacet[j].Points[n];
												double z_Q = co_pTool[Bearing].pPCoord[N][2];
												if (z_Q < zM) {
													zM = z_Q;
													M = co_pTool[Bearing].pFacet[j].Points[n];        // Find node with lower z //
												}
											}
											for (int e = j; e < co_pTool[Bearing].nFacets; e++) {       // remove element   //
												if (e != co_pTool[Bearing].nFacets) {
													co_pTool[Bearing].pFacet[e].Points[0] = co_pTool[Bearing].pFacet[e + 1].Points[0];
													co_pTool[Bearing].pFacet[e].Points[1] = co_pTool[Bearing].pFacet[e + 1].Points[1];
													co_pTool[Bearing].pFacet[e].Points[2] = co_pTool[Bearing].pFacet[e + 1].Points[2];
												}
											}

											if ((co_pTool[Bearing].pFacet = (CO_FACET_DATA *)GE_REALLOC(CO_FACET_DATA, [1],
												co_pTool[Bearing].pFacet, co_pTool[Bearing].nFacets - 1)) == NULL)
												return 0;

											co_pTool[Bearing].nFacets -= 1;
											j = 0;
											k = 0;

											//if (!rm_co_ToolInitialisation(12, co_pTool + Bearing)) return 0;
											rm_co_NodeDataInitialisation(me_nBNodes, co_pBNdata);

											for (int l = 0; l < co_pTool[Bearing].nFacets; l++) {      // remove elements that have the node with lower z   //
												for (int q = 0; q < 3; q++) {
													int Q = co_pTool[Bearing].pFacet[l].Points[q];
													if (Q == M) {
														int remove = 0;                                       // avoid removing last elements  //
														for (int e = 0; e < 3; e++) {
															int E = co_pTool[Bearing].pFacet[l].Points[q];
															double z_E = co_pTool[Bearing].pPCoord[E][2];
															if (z_E > z_min+0.1) remove++;
														}
														if (remove == 3) {
															for (int p = l; p < co_pTool[Bearing].nFacets; p++) {       // remove element   //
																if (p != co_pTool[Bearing].nFacets) {
																	co_pTool[Bearing].pFacet[p].Points[0] = co_pTool[Bearing].pFacet[p + 1].Points[0];
																	co_pTool[Bearing].pFacet[p].Points[1] = co_pTool[Bearing].pFacet[p + 1].Points[1];
																	co_pTool[Bearing].pFacet[p].Points[2] = co_pTool[Bearing].pFacet[p + 1].Points[2];
																}
															}

															if ((co_pTool[Bearing].pFacet = (CO_FACET_DATA *)GE_REALLOC(CO_FACET_DATA, [1],
																co_pTool[Bearing].pFacet, co_pTool[Bearing].nFacets - 1)) == NULL)
																return 0;

															co_pTool[Bearing].nFacets -= 1;
															l = 0;
															q = 0;

															//if (!rm_co_ToolInitialisation(12, co_pTool + Bearing)) return 0;
															rm_co_NodeDataInitialisation(me_nBNodes, co_pBNdata);
														}
													}
												}
											}
										}
									}
								}


								// UPDATE me_pTriaReg
								for (int m = 0; m < me_pTriaReg[reg].nTrias; m++) {
									int remove = 0;
									for (int e = 0; e < 3; e++) {
										double z_A = me_pNdata[me_pTriaReg[reg].pTriaNodes[m][e]].OldCoord[2];
										if (z_A > z_min+0.1) remove++;
									}

									if (remove == 3) {      // do not remove last elements
										int F = me_pTriaReg[reg].pTriaNodes[m][0];
										double zF = me_pNdata[me_pTriaReg[reg].pTriaNodes[m][0]].OldCoord[2];
										for (int t = 0; t < 3; t++) {
											double zD = me_pNdata[me_pTriaReg[reg].pTriaNodes[m][t]].OldCoord[2];
											if (zD < zF) {
												F = me_pTriaReg[reg].pTriaNodes[m][t];
												zF = zD;
											}
										}
										for (int n = 0; n < 3; n++) {
											//remove elements containing node
											if (me_pTriaReg[reg].pTriaNodes[m][n] == BearingOptInfo[i].ID) {
												for (int f = m; f < me_pTriaReg[reg].nTrias; f++) {
													if (f != me_pTriaReg[reg].nTrias) {
														me_pTriaReg[reg].pTriaNodes[f][0] = me_pTriaReg[reg].pTriaNodes[f + 1][0];
														me_pTriaReg[reg].pTriaNodes[f][1] = me_pTriaReg[reg].pTriaNodes[f + 1][1];
														me_pTriaReg[reg].pTriaNodes[f][2] = me_pTriaReg[reg].pTriaNodes[f + 1][2];
													}
												}

												if ((me_pTriaReg[reg].pTriaNodes = (int(*)[3])GE_REALLOC(int, [3],
													me_pTriaReg[reg].pTriaNodes, me_pTriaReg[reg].nTrias - 1)) == NULL)
													return 0;

												if ((me_pTriaReg[reg].BearingLine = (int(*)[3])GE_REALLOC(int, [3],
													me_pTriaReg[reg].BearingLine, me_pTriaReg[reg].nTrias - 1)) == NULL)
													return 0;

												if ((me_pTriaReg[reg].BearingElement = (int(*)[1])GE_REALLOC(int, [1],
													me_pTriaReg[reg].BearingElement, me_pTriaReg[reg].nTrias - 1)) == NULL)
													return 0;

												me_pTriaReg[reg].nTrias -= 1;
												m = 0;
												n = 0;

												for (int w = 0; w < me_pTriaReg[reg].nTrias; w++) {         // remove elements that have the node with lower z   //
													int remove = 0;
													for (int g = 0; g < 3; g++) {
														double zG = me_pNdata[me_pTriaReg[reg].pTriaNodes[w][g]].OldCoord[2];
														if (zG > z_min+0.1) remove++;
													}
													if (remove == 3) {      // avoid removing last element  //
														for (int s = 0; s < 3; s++) {
															if (me_pTriaReg[reg].pTriaNodes[w][s] == F) {
																for (int f = w; f < me_pTriaReg[reg].nTrias; f++) {
																	if (f != me_pTriaReg[reg].nTrias) {
																		me_pTriaReg[reg].pTriaNodes[f][0] = me_pTriaReg[reg].pTriaNodes[f + 1][0];
																		me_pTriaReg[reg].pTriaNodes[f][1] = me_pTriaReg[reg].pTriaNodes[f + 1][1];
																		me_pTriaReg[reg].pTriaNodes[f][2] = me_pTriaReg[reg].pTriaNodes[f + 1][2];
																	}
																}

																if ((me_pTriaReg[reg].pTriaNodes = (int(*)[3])GE_REALLOC(int, [3],
																	me_pTriaReg[reg].pTriaNodes, me_pTriaReg[reg].nTrias - 1)) == NULL)
																	return 0;

																if ((me_pTriaReg[reg].BearingLine = (int(*)[3])GE_REALLOC(int, [3],
																	me_pTriaReg[reg].BearingLine, me_pTriaReg[reg].nTrias - 1)) == NULL)
																	return 0;

																if ((me_pTriaReg[reg].BearingElement = (int(*)[1])GE_REALLOC(int, [1],
																	me_pTriaReg[reg].BearingElement, me_pTriaReg[reg].nTrias - 1)) == NULL)
																	return 0;

																me_pTriaReg[reg].nTrias -= 1;
																w = 0;
																s = 0;
															}
														}
													}
												}
											}
										}

									}
								}

							}
						}
					}
				}
			}
		}
	}
					


	/*--------------------------------------------------------------------------------------------------------*/
	//                            STRETCHING BEARING LINE   VRD > 1
	/*--------------------------------------------------------------------------------------------------------*/
	
	if (debug == 1) {
		FILE *Increase_Nodes;

		char vtkFileName[255];

		sprintf(vtkFileName, "Increase Nodes_%i.csv", nl_counterAndrea);

		Increase_Nodes = fopen(vtkFileName, "w");
		fprintf(Increase_Nodes, "x,y,z,\n");
		for (int i = 0; i < nFoundNodes; i++) {
			if (BearingOptInfo[i].BearNodesVDR > 1+tol) {
				double x = me_pNdata[BearingOptInfo[i].ID].OldCoord[0];
				double y = me_pNdata[BearingOptInfo[i].ID].OldCoord[1];
				double z = me_pNdata[BearingOptInfo[i].ID].OldCoord[2];
				fprintf(Increase_Nodes, "%f,%f,%f\n", x, y, z);
			}
		}
		fclose(Increase_Nodes);
	}


	if (regionID > -1) {
	for (int reg = 0; reg < me_nTriaReg; reg++) {
		if (nl_pPdata[0].pBearingRegionFlag[reg] == 1) {
			int Bearing = -1;
			for (int m = 0; m < co_nTools; m++) {
				if (strcmp(me_pTriaReg[reg].Name, co_pTool[m].ToolName) == 0) {
					Bearing = m;
				}
			}
			// add Trias with two bearing nodes
			for (int i = 0; i < me_nStriangles; i++) {
				int added = 0;
				double x[3];
				double y[3];
				double z[3];
				int el[3];
				el[0] = -1;
				el[1] = -1;
				el[2] = -1;
				int found = 0;
				for (int k = 0; k < 3; k++) {
					x[k] = me_pNdata[me_pSTdata[i].Nodes[k]].OldCoord[0];
					y[k] = me_pNdata[me_pSTdata[i].Nodes[k]].OldCoord[1];
					z[k] = me_pNdata[me_pSTdata[i].Nodes[k]].OldCoord[2];



					if (z[k] >= z_min+100) found = found + 10;                // limit max bearing length

					for (int j = 0; j < nFoundNodes; j++) {
						if (BearingOptInfo[j].ID == me_pSTdata[i].Nodes[k] && BearingOptInfo[j].BearNodesVDR > 1+tol) {
							el[found] = j;
							found++;
						}
					}
				}

				if (found == 2) {
					int add = 1;
					for (int m = 0; m < co_pTool[Bearing].nFacets; m++) {
						int ee = 0;
						for (int n = 0; n < 3; n++) {
							int P = co_pTool[Bearing].pFacet[m].Points[n];
							double x_P = co_pTool[Bearing].pPCoord[P][0];
							double y_P = co_pTool[Bearing].pPCoord[P][1];
							double z_P = co_pTool[Bearing].pPCoord[P][2];
							for (int k = 0; k < 3; k++) {
								if (fabs(x_P - x[k]) < tol2 && fabs(y_P - y[k]) < tol2 && fabs(z_P - z[k]) < tol2) ee++;
							}
						}
						if (ee == 3) add = 0;
					}
					if (add == 1) {
						int nod[3];
						nod[0] = -1;
						nod[1] = -1;
						nod[2] = -1;
						int ele[3];
						ele[0] = -1;
						ele[1] = -1;
						ele[2] = -1;
						int share = 0;
						for (int m = 0; m < co_pTool[Bearing].nFacets; m++) {
							double x_P[3];
							double y_P[3];
							double z_P[3];
							for (int n = 0; n < 3; n++) {
								int P = co_pTool[Bearing].pFacet[m].Points[n];
								x_P[n] = co_pTool[Bearing].pPCoord[P][0];
								y_P[n] = co_pTool[Bearing].pPCoord[P][1];
								z_P[n] = co_pTool[Bearing].pPCoord[P][2];
								for (int k = 0; k < 3; k++) {
									if (fabs(x_P[n] - x[k]) < tol2 && fabs(y_P[n] - y[k]) < tol2 && fabs(z_P[n] - z[k]) < tol2
										&& nod[k] == -1) {
										nod[k] = n;
										ele[k] = m;
										share++;
									}
								}
							}
						}
						if (share == 2 && added == 0) {
							
							if ((co_pTool[Bearing].pFacet = (CO_FACET_DATA *)GE_REALLOC(CO_FACET_DATA, [1],
								co_pTool[Bearing].pFacet, co_pTool[Bearing].nFacets + 1)) == NULL)
								return 0;

							if ((co_pTool[Bearing].pPCoord = (double(*)[3])GE_REALLOC(double, [3],
								co_pTool[Bearing].pPCoord, co_pTool[Bearing].nPoints + 1)) == NULL)
								return 0;

							// UPDATE me_pTriaReg
							if ((me_pTriaReg[reg].pTriaNodes = (int(*)[3])GE_REALLOC(int, [3],
								me_pTriaReg[reg].pTriaNodes, me_pTriaReg[reg].nTrias + 1)) == NULL)
								return 0;

							if ((me_pTriaReg[reg].pNodes = (int(*))GE_REALLOC(int, [1],
								me_pTriaReg[reg].pNodes, me_pTriaReg[reg].nNodes + 1)) == NULL)
								return 0;

							if ((me_pTriaReg[reg].BearingLine = (int(*)[3])GE_REALLOC(int, [3],
								me_pTriaReg[reg].BearingLine, me_pTriaReg[reg].nTrias + 1)) == NULL)
								return 0;

							if ((me_pTriaReg[reg].BearingElement = (int(*)[1])GE_REALLOC(int, [1],
								me_pTriaReg[reg].BearingElement, me_pTriaReg[reg].nTrias + 1)) == NULL)
								return 0;
							

							
							int b = 0;
							for (int a = 0; a < 3; a++) {
								if (nod[a] == -1) {
									co_pTool[Bearing].pPCoord[co_pTool[Bearing].nPoints][0] = x[a];
									co_pTool[Bearing].pPCoord[co_pTool[Bearing].nPoints][1] = y[a];
									co_pTool[Bearing].pPCoord[co_pTool[Bearing].nPoints][2] = z[a];
									co_pTool[Bearing].pFacet[co_pTool[Bearing].nFacets].Points[2] = co_pTool[Bearing].nPoints;

									// UPDATE me_pTriaReg 
									me_pTriaReg[reg].pNodes[me_pTriaReg[reg].nNodes] = me_pSTdata[i].Nodes[a];
									me_pTriaReg[reg].pTriaNodes[me_pTriaReg[reg].nTrias][2] = me_pSTdata[i].Nodes[a];

								}
								if (nod[a] != -1) {
									co_pTool[Bearing].pFacet[co_pTool[Bearing].nFacets].Points[b] = co_pTool[Bearing].pFacet[ele[a]].Points[nod[a]];

									// UPDATE me_pTriaReg 
									me_pTriaReg[reg].pTriaNodes[me_pTriaReg[reg].nTrias][b] = me_pSTdata[i].Nodes[a];
									b++;
								}
							}

							added = 1;
							co_pTool[Bearing].nFacets += 1;
							co_pTool[Bearing].nPoints += 1;

							// UPDATE me_pTriaReg 
							me_pTriaReg[reg].nTrias += 1;
							me_pTriaReg[reg].nNodes += 1;

							terminate = 0;
						}
						if (share == 3 && added == 0) {
							
							if ((co_pTool[Bearing].pFacet = (CO_FACET_DATA *)GE_REALLOC(CO_FACET_DATA, [1],
								co_pTool[Bearing].pFacet, co_pTool[Bearing].nFacets + 1)) == NULL)
								return 0;

							// UPDATE me_pTriaReg
							if ((me_pTriaReg[reg].pTriaNodes = (int(*)[3])GE_REALLOC(int, [3],
								me_pTriaReg[reg].pTriaNodes, me_pTriaReg[reg].nTrias + 1)) == NULL)
								return 0;


							if ((me_pTriaReg[reg].BearingLine = (int(*)[3])GE_REALLOC(int, [3],
								me_pTriaReg[reg].BearingLine, me_pTriaReg[reg].nTrias + 1)) == NULL)
								return 0;

							if ((me_pTriaReg[reg].BearingElement = (int(*)[1])GE_REALLOC(int, [1],
								me_pTriaReg[reg].BearingElement, me_pTriaReg[reg].nTrias + 1)) == NULL)
								return 0;
						

							for (int a = 0; a < 3; a++) {

								co_pTool[Bearing].pFacet[co_pTool[Bearing].nFacets].Points[a] = co_pTool[Bearing].pFacet[ele[a]].Points[nod[a]];

								// UPDATE me_pTriaReg 
								me_pTriaReg[reg].pTriaNodes[me_pTriaReg[reg].nTrias][a] = me_pSTdata[i].Nodes[a];

							}

							added = 1;
							co_pTool[Bearing].nFacets += 1;
							terminate = 0;

							// UPDATE me_pTriaReg 
							me_pTriaReg[reg].nTrias += 1;

						}
					
					}
				}

			}
			
			// add trias with one bearing node 
			for (int i = 0; i < me_nStriangles; i++) {
				int added = 0;
				double x[3];
				double y[3];
				double z[3];
				int e[3];
				e[0] = -1;
				e[1] = -1;
				e[2] = -1;
				int found = 0;
				for (int k = 0; k < 3; k++) {
					x[k] = me_pNdata[me_pSTdata[i].Nodes[k]].OldCoord[0];
					y[k] = me_pNdata[me_pSTdata[i].Nodes[k]].OldCoord[1];
					z[k] = me_pNdata[me_pSTdata[i].Nodes[k]].OldCoord[2];

					if (z[k] >= 100) found = found + 10;                // limit max bearing length

					for (int j = 0; j < nFoundNodes; j++) {
						if (BearingOptInfo[j].ID == me_pSTdata[i].Nodes[k] && BearingOptInfo[j].BearNodesVDR > 1+tol) {
							e[found] = j;
							found++;
						}
					}
				}
				for (int k = 0; k < 3; k++) {
					for (int l = 0; l < nFoundNodes; l++) {
						if (BearingOptInfo[l].BearNodesVDR < 1+tol && x[k] == BearingOptInfo[l].BearNodesXYZ[0] &&
							y[k] == BearingOptInfo[l].BearNodesXYZ[1] && 
							BearingOptInfo[l].ClosestBearingNode == 1) found = found + 10;
					}
				}
				if (found == 1) {
					int add = 1;
					for (int m = 0; m < co_pTool[Bearing].nFacets; m++) {
						int ee = 0;
						for (int n = 0; n < 3; n++) {
							int P = co_pTool[Bearing].pFacet[m].Points[n];
							double x_P = co_pTool[Bearing].pPCoord[P][0];
							double y_P = co_pTool[Bearing].pPCoord[P][1];
							double z_P = co_pTool[Bearing].pPCoord[P][2];
							for (int k = 0; k < 3; k++) {
								if (fabs(x_P - x[k]) < tol2 && fabs(y_P - y[k]) < tol2 && fabs(z_P - z[k]) < tol2) ee++;
							}
						}
						if (ee == 3) add = 0;
					}
					if (add == 1) {
						int share = 0;
						int nod[3];
						nod[0] = -1;
						nod[1] = -1;
						nod[2] = -1;
						int el[3];
						el[0] = -1;
						el[1] = -1;
						el[2] = -1;
						for (int m = 0; m < co_pTool[Bearing].nFacets; m++) {
							double x_P[3];
							double y_P[3];
							double z_P[3];
							for (int n = 0; n < 3; n++) {
								int P = co_pTool[Bearing].pFacet[m].Points[n];
								x_P[n] = co_pTool[Bearing].pPCoord[P][0];
								y_P[n] = co_pTool[Bearing].pPCoord[P][1];
								z_P[n] = co_pTool[Bearing].pPCoord[P][2];
								for (int k = 0; k < 3; k++) {
									if (fabs(x_P[n] - x[k]) < tol2 && fabs(y_P[n] - y[k]) < tol2 && fabs(z_P[n] - z[k]) < tol2
										&& nod[k] == -1) {
										nod[k] = n;
										el[k] = m;
										share++;
									}
								}
							}
						}
						if (share == 2 && added == 0) {
						
							if ((co_pTool[Bearing].pFacet = (CO_FACET_DATA *)GE_REALLOC(CO_FACET_DATA, [1],
								co_pTool[Bearing].pFacet, co_pTool[Bearing].nFacets + 1)) == NULL)
								return 0;

							if ((co_pTool[Bearing].pPCoord = (double(*)[3])GE_REALLOC(double, [3],
								co_pTool[Bearing].pPCoord, co_pTool[Bearing].nPoints + 1)) == NULL)
								return 0;

							// UPDATE me_pTriaReg 
							if ((me_pTriaReg[reg].pTriaNodes = (int(*)[3])GE_REALLOC(int, [3],
								me_pTriaReg[reg].pTriaNodes, me_pTriaReg[reg].nTrias + 1)) == NULL)
								return 0;

							if ((me_pTriaReg[reg].pNodes = (int(*))GE_REALLOC(int, [1],
								me_pTriaReg[reg].pNodes, me_pTriaReg[reg].nNodes + 1)) == NULL)
								return 0;

							if ((me_pTriaReg[reg].BearingLine = (int(*)[3])GE_REALLOC(int, [3],
								me_pTriaReg[reg].BearingLine, me_pTriaReg[reg].nTrias + 1)) == NULL)
								return 0;

							if ((me_pTriaReg[reg].BearingElement = (int(*)[1])GE_REALLOC(int, [1],
								me_pTriaReg[reg].BearingElement, me_pTriaReg[reg].nTrias + 1)) == NULL)
								return 0;

							


							int b = 0;
							for (int a = 0; a < 3; a++) {
								if (nod[a] == -1) {
									co_pTool[Bearing].pPCoord[co_pTool[Bearing].nPoints][0] = x[a];
									co_pTool[Bearing].pPCoord[co_pTool[Bearing].nPoints][1] = y[a];
									co_pTool[Bearing].pPCoord[co_pTool[Bearing].nPoints][2] = z[a];
									co_pTool[Bearing].pFacet[co_pTool[Bearing].nFacets].Points[2] = co_pTool[Bearing].nPoints;

									// UPDATE me_pTriaReg 
									me_pTriaReg[reg].pNodes[me_pTriaReg[reg].nNodes] = me_pSTdata[i].Nodes[a];
									me_pTriaReg[reg].pTriaNodes[me_pTriaReg[reg].nTrias][2] = me_pSTdata[i].Nodes[a];
								}
								if (nod[a] != -1) {
									co_pTool[Bearing].pFacet[co_pTool[Bearing].nFacets].Points[b] = co_pTool[Bearing].pFacet[el[a]].Points[nod[a]];

									// UPDATE me_pTriaReg 
									me_pTriaReg[reg].pTriaNodes[me_pTriaReg[reg].nTrias][b] = me_pSTdata[i].Nodes[a];
									b++;
								}
							}

							added = 1;
							co_pTool[Bearing].nFacets += 1;
							co_pTool[Bearing].nPoints += 1;

							// UPDATE me_pTriaReg 
							me_pTriaReg[reg].nTrias += 1;
							me_pTriaReg[reg].nNodes += 1;

							terminate = 0;

						}

						if (share == 3 && added == 0) {

							if ((co_pTool[Bearing].pFacet = (CO_FACET_DATA *)GE_REALLOC(CO_FACET_DATA, [1],
								co_pTool[Bearing].pFacet, co_pTool[Bearing].nFacets + 1)) == NULL)
								return 0;

							// UPDATE me_pTriaReg 
							if ((me_pTriaReg[reg].pTriaNodes = (int(*)[3])GE_REALLOC(int, [3],
								me_pTriaReg[reg].pTriaNodes, me_pTriaReg[reg].nTrias + 1)) == NULL)
								return 0;

							if ((me_pTriaReg[reg].BearingLine = (int(*)[3])GE_REALLOC(int, [3],
								me_pTriaReg[reg].BearingLine, me_pTriaReg[reg].nTrias + 1)) == NULL)
								return 0;

							if ((me_pTriaReg[reg].BearingElement = (int(*)[1])GE_REALLOC(int, [1],
								me_pTriaReg[reg].BearingElement, me_pTriaReg[reg].nTrias + 1)) == NULL)
								return 0;

							for (int a = 0; a < 3; a++) {
								co_pTool[Bearing].pFacet[co_pTool[Bearing].nFacets].Points[a] = co_pTool[Bearing].pFacet[el[a]].Points[nod[a]];

								// UPDATE me_pTriaReg 
								me_pTriaReg[reg].pTriaNodes[me_pTriaReg[reg].nTrias][a] = me_pSTdata[i].Nodes[a];
							}

							added = 1;
							co_pTool[Bearing].nFacets += 1;

							// UPDATE me_pTriaReg 
							me_pTriaReg[reg].nTrias += 1;

							terminate = 0;
						}
					}
				}
			}
				
			/*--------------------------------------------------------------------------------------------------------*/
			//                                              VOID FILLING
			/*--------------------------------------------------------------------------------------------------------*/
			for (int i = 0; i < me_nStriangles; i++) {
				int added = 0;
				int add = 1;
				double x[3];
				double y[3];
				double z[3];
				int found = 0;
				for (int k = 0; k < 3; k++) {
					x[k] = me_pNdata[me_pSTdata[i].Nodes[k]].OldCoord[0];
					y[k] = me_pNdata[me_pSTdata[i].Nodes[k]].OldCoord[1];
					z[k] = me_pNdata[me_pSTdata[i].Nodes[k]].OldCoord[2];
				}
				for (int m = 0; m < co_pTool[Bearing].nFacets; m++) {
					int ee = 0;
					for (int n = 0; n < 3; n++) {
						int P = co_pTool[Bearing].pFacet[m].Points[n];
						double x_P = co_pTool[Bearing].pPCoord[P][0];
						double y_P = co_pTool[Bearing].pPCoord[P][1];
						double z_P = co_pTool[Bearing].pPCoord[P][2];
						for (int k = 0; k < 3; k++) {
							if (fabs(x_P - x[k]) < tol2 && fabs(y_P - y[k]) < tol2 && fabs(z_P - z[k]) < tol2) ee++;
						}
					}
					if (ee == 3) add = 0;
				}
				if (add == 1) {
					int share = 0;
					int nod[3];
					nod[0] = -1;
					nod[1] = -1;
					nod[2] = -1;
					int el[3];
					el[0] = -1;
					el[1] = -1;
					el[2] = -1;
					for (int m = 0; m < co_pTool[Bearing].nFacets; m++) {
						double x_P[3];
						double y_P[3];
						double z_P[3];
						for (int n = 0; n < 3; n++) {
							int P = co_pTool[Bearing].pFacet[m].Points[n];
							x_P[n] = co_pTool[Bearing].pPCoord[P][0];
							y_P[n] = co_pTool[Bearing].pPCoord[P][1];
							z_P[n] = co_pTool[Bearing].pPCoord[P][2];
							for (int k = 0; k < 3; k++) {
								if (fabs(x_P[n] - x[k]) < tol2 && fabs(y_P[n] - y[k]) < tol2 && fabs(z_P[n] - z[k]) < tol2
									&& nod[k] == -1) {
									nod[k] = n;
									el[k] = m;
									share++;
								}
							}
						}
					}
					if (share == 3 && added == 0) {

						if ((co_pTool[Bearing].pFacet = (CO_FACET_DATA *)GE_REALLOC(CO_FACET_DATA, [1],
							co_pTool[Bearing].pFacet, co_pTool[Bearing].nFacets + 1)) == NULL)
							return 0;

						// UPDATE me_pTriaReg 
						if ((me_pTriaReg[reg].pTriaNodes = (int(*)[3])GE_REALLOC(int, [3],
							me_pTriaReg[reg].pTriaNodes, me_pTriaReg[reg].nTrias + 1)) == NULL)
							return 0;

						if ((me_pTriaReg[reg].BearingLine = (int(*)[3])GE_REALLOC(int, [3],
							me_pTriaReg[reg].BearingLine, me_pTriaReg[reg].nTrias + 1)) == NULL)
							return 0;

						if ((me_pTriaReg[reg].BearingElement = (int(*)[1])GE_REALLOC(int, [1],
							me_pTriaReg[reg].BearingElement, me_pTriaReg[reg].nTrias + 1)) == NULL)
							return 0;

						for (int a = 0; a < 3; a++) {
							co_pTool[Bearing].pFacet[co_pTool[Bearing].nFacets].Points[a] = co_pTool[Bearing].pFacet[el[a]].Points[nod[a]];

							// UPDATE me_pTriaReg 
							me_pTriaReg[reg].pTriaNodes[me_pTriaReg[reg].nTrias][a] = me_pSTdata[i].Nodes[a];
						}

						added = 1;
						co_pTool[Bearing].nFacets += 1;

						// UPDATE me_pTriaReg 
						me_pTriaReg[reg].nTrias += 1;

					}
				}
			}

			/*--------------------------------------------------------------------------------------------------------*/
			//                            TOOL RE-BUILT AND VARIABLES RE-INITIALIZATION
			/*--------------------------------------------------------------------------------------------------------*/

			co_pTool[Bearing].nEdges = 0;


			if (debug == 1) {
				char vtkFileName[255];
				FILE *vtkFile;



				sprintf(vtkFileName, "%s_Optimized_%i.vtk", co_pTool[Bearing].ToolName, nl_counterAndrea);

				vtkFile = fopen(vtkFileName, "w");
				fprintf(vtkFile, "# vtk DataFile Version 2.0");
				fprintf(vtkFile, "\nExtrusion Processtime %f", nl_counterAndrea);
				fprintf(vtkFile, "\nASCII");

				fprintf(vtkFile, "\n\nDATASET UNSTRUCTURED_GRID");
				fprintf(vtkFile, "\nPOINTS %d float", co_pTool[Bearing].nPoints);

				for (int no = 0; no < co_pTool[Bearing].nPoints; no++)
				{
					fprintf(vtkFile, "\n%f %f %f", co_pTool[Bearing].pPCoord[no][0],
						co_pTool[Bearing].pPCoord[no][1],
						co_pTool[Bearing].pPCoord[no][2]);
				}

				fprintf(vtkFile, "\n\nCELLS %d %d", co_pTool[Bearing].nFacets, 4 * co_pTool[Bearing].nFacets);
				for (int el = 0; el < co_pTool[Bearing].nFacets; el++)
				{
					fprintf(vtkFile, "\n%d %d %d %d", 3, co_pTool[Bearing].pFacet[el].Points[0],
						co_pTool[Bearing].pFacet[el].Points[1],
						co_pTool[Bearing].pFacet[el].Points[2]);
				}

				fprintf(vtkFile, "\n\nCELL_TYPES %d", co_pTool[Bearing].nFacets);
				for (int il = 0; il < co_pTool[Bearing].nFacets; il++)
				{
					fprintf(vtkFile, "\n5");
				}

				fprintf(vtkFile, "\n");

				fclose(vtkFile);
			}

			if (debug == 1) {
				char vtkFileName[255];
				FILE *vtkFile;


				sprintf(vtkFileName, "%s_Optimized_me_pTria_%i.vtk", co_pTool[Bearing].ToolName, nl_counterAndrea);

				vtkFile = fopen(vtkFileName, "w");
				fprintf(vtkFile, "# vtk DataFile Version 2.0");
				fprintf(vtkFile, "\nExtrusion Processtime %f", counter);
				fprintf(vtkFile, "\nASCII");

				fprintf(vtkFile, "\n\nDATASET UNSTRUCTURED_GRID");
				fprintf(vtkFile, "\nPOINTS %d float", me_pTriaReg[regionID].nNodes);

				for (int no = 0; no < me_pTriaReg[regionID].nNodes; no++) {

					fprintf(vtkFile, "\n%f %f %f", me_pNdata[me_pTriaReg[regionID].pNodes[no]].OldCoord[0],
						me_pNdata[me_pTriaReg[regionID].pNodes[no]].OldCoord[1],
						me_pNdata[me_pTriaReg[regionID].pNodes[no]].OldCoord[2]);

				}

				fprintf(vtkFile, "\n\nCELLS %d %d", me_pTriaReg[regionID].nTrias, 4 * me_pTriaReg[regionID].nTrias);
				for (int el = 0; el < me_pTriaReg[regionID].nTrias; el++) {
					int a[3];
					int q = 0;
					for (int k = 0; k < 3; k++) {
						int added = 0;
						for (int no = 0; no < me_pTriaReg[regionID].nNodes; no++) {
							if (added == 0 
								&& me_pTriaReg[regionID].pTriaNodes[el][k] == me_pTriaReg[regionID].pNodes[no]) {

								a[q] = no;
								added = 1;
								q++;

							}
						}
					}
					fprintf(vtkFile, "\n%d %d %d %d", 3, a[0],
						a[1],
						a[2]);
				}

				fprintf(vtkFile, "\n\nCELL_TYPES %d", me_pTriaReg[regionID].nTrias);
				for (int il = 0; il < me_pTriaReg[regionID].nTrias; il++)
				{
					fprintf(vtkFile, "\n5");
				}

				fprintf(vtkFile, "\n");

				fclose(vtkFile);
			}



			/*if (!rm_co_ToolInitialisation(12, co_pTool + Bearing)) return 0;
			rm_co_NodeDataInitialisation(me_nBNodes, co_pBNdata);*/

			for (int i = 0; i < me_pTriaReg[reg].nTrias; i++)
			{
				me_pTriaReg[reg].BearingLine[i][0] = -1;
				me_pTriaReg[reg].BearingLine[i][1] = -1;
				me_pTriaReg[reg].BearingLine[i][2] = -1;
				me_pTriaReg[reg].BearingElement[i] = -1;
			}
			
		}
	}
	
}



	for (int j = 0; j < rm_nTTrias; j++) {
		for (int k = 0; k < 3; k++) {
			rm_pTTdata[j].checked[k] = 0;
			rm_pTTdata[j].ClosestBearingNode[k] = -1;
		}
	}



/*--------------------------------------------------------------------------------------------------------*/
//                             PRINT UPDATED MESH && TERMINATING CONDITION
/*--------------------------------------------------------------------------------------------------------*/

	char pfFileName2[255];
	FILE *pfFile2;
	sprintf(pfFileName2, "mesh_Optimized.pf");

	pfFile2 = fopen(pfFileName2, "w");
	// Write new mesh
	fprintf(pfFile2, "Nodes\n");

	for (int i = 0; i < me_nNodes; i++)
	{
		fprintf(pfFile2, "%i %.6f %.6f %.6f \n", i, me_pNdata[i].OldCoord[0], me_pNdata[i].OldCoord[1], me_pNdata[i].OldCoord[2]);
	}

	fprintf(pfFile2, "\nTetra\n");

	for (int i = 0; i < me_nTetras; i++)
	{
		fprintf(pfFile2, "%i %i %i %i %i \n", i, me_pETdata[i].Nodes[0], me_pETdata[i].Nodes[1], me_pETdata[i].Nodes[2], me_pETdata[i].Nodes[3]);
	}

	for (int m = 0; m < me_nTriaReg; m++) {
		 fprintf(pfFile2, "\nFacets %s \n", me_pTriaReg[m].Name);
		 for (int i = 0; i < me_pTriaReg[m].nTrias;i++)
		 {
			 fprintf(pfFile2, "%i %i %i %i \n", i, me_pTriaReg[m].pTriaNodes[i][0], me_pTriaReg[m].pTriaNodes[i][1], me_pTriaReg[m].pTriaNodes[i][2]);
		 }
	}

	fclose(pfFile2);


	if (terminate == 1) {
		char pfFileName2[255];
		FILE *term;
		sprintf(term, "terminate.txt");
		fclose(term);
	}
	free(BearingOptInfo);
	return 0;
}