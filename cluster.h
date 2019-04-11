#include"octree.h"
#define POS_SCALE 5
#define NOR_SCALE 1
#define CLUSTER_SIZE 32

float distance_pos(struct Intersection sp1, struct Intersection sp2, float scale)//位置
{

	float res, distance_nor;
	res = 0;
	res = sqrt((sp1.p[0] - sp2.p[0])*(sp1.p[0] - sp2.p[0]) + (sp1.p[1] - sp2.p[1])*(sp1.p[1] - sp2.p[1]) + (sp1.p[2] - sp2.p[2])*(sp1.p[2] - sp2.p[2]));
	return res;
}
float distance_nor(struct Intersection sp1, struct Intersection sp2, float scale)//法向
{
	float res;
	res = sqrt((sp1.n[0] - sp2.n[0])*(sp1.n[0] - sp2.n[0]) + (sp1.n[1] - sp2.n[1])*(sp1.n[1] - sp2.n[1]) + (sp1.n[2] - sp2.n[2])*(sp1.n[2] - sp2.n[2]));
	res *= scale;
	res += distance_pos(sp1, sp2, scale);
	return res;
}
void swap_sp(struct Intersection *sps, int i, int j)
{
	struct Intersection sp_temp;
	sp_temp = sps[i];
	sps[i] = sps[j];
	sps[j] = sp_temp;
}
void swap_index(int *index_list, int i, int j){
	int temp;
	temp = index_list[i];
	index_list[i] = index_list[j];
	index_list[j] = temp;
}
void center(struct Intersection *sps, struct Intersection *sp0, int i, int n_16)
{
	sp0->p[0] = 0;
	sp0->p[1] = 0;
	sp0->p[2] = 0;
	sp0->n[0] = 0;
	sp0->n[1] = 0;
	sp0->n[2] = 0;
	int j;
	for (j = i; j<i + n_16; j++)
	{
		sp0->p[0] += sps[j].p[0];
		sp0->p[1] += sps[j].p[1];
		sp0->p[2] += sps[j].p[2];
		sp0->n[0] += sps[j].n[0];
		sp0->n[1] += sps[j].n[1];
		sp0->n[2] += sps[j].n[2];
	}
	sp0->p[0] = sp0->p[0] / n_16;
	sp0->p[1] = sp0->p[1] / n_16;
	sp0->p[2] = sp0->p[2] / n_16;
	sp0->n[0] = sp0->n[0] / n_16;
	sp0->n[1] = sp0->n[1] / n_16;
	sp0->n[2] = sp0->n[2] / n_16;

}

void variance(struct Intersection *sps, struct Intersection sp0, int i, float var[CLUSTER_SIZE], float *maxvar, int *varnum, float scale)//计算16个元素方差求出最大的一�?
{
	*varnum = 0;
	*maxvar = 0;
	int j;
	for ( j = i; j<i + CLUSTER_SIZE; j++)
	{
		var[j - i] = distance_nor(sp0, sps[j], scale);
		if (var[j - i]>*maxvar)
		{
			*maxvar = var[j - i];
			*varnum = j;
		}
	}
}

void cluster(struct Intersection *sps, int *task_index, int num_ray, float scale) 
{
	struct Intersection sp0; //质心
	float dis_pos, dis_nor;
	int n = num_ray / CLUSTER_SIZE;
	int first, second, n_16;
	float var[CLUSTER_SIZE];  //一组元素的方差
	float maxvar; //最大方�?
	int varnum;    //最大方差所在的位置
	int changed = 0;
	int i,j;
	for (i = 1; i<n; i++)
	{
		changed = 0;
		center(sps, &sp0, (i - 1)*CLUSTER_SIZE, CLUSTER_SIZE);
		variance(sps, sp0, (i - 1)*CLUSTER_SIZE, var, &maxvar, &varnum, scale);
		//int nrays = ((i*16+300)>num_ray)?num_ray:(i*16+200);
		for ( j = i*CLUSTER_SIZE; j<num_ray; j++)
		{
			if (changed >= 300)
				continue;
			float new_nor = distance_nor(sp0, sps[j], scale);
			float new_pos = distance_pos(sp0, sps[j], scale);
			if (new_nor<maxvar){
				swap_sp(sps, varnum, j);
				swap_index(task_index, varnum, j);
				center(sps, &sp0, (i - 1)*CLUSTER_SIZE, CLUSTER_SIZE);
				variance(sps, sp0, (i - 1)*CLUSTER_SIZE, var, &maxvar, &varnum, scale);
				changed = 0;
			}
			changed++;
		}
	}
}
