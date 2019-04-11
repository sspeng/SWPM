#include <stdio.h>
#include <math.h>
#include <string.h>
#include "slave.h"

#define M_PI 3.14159265358979323846
#define INV_PI 0.31830988618379067154
#define inline __inline
#define PACKET_SIZE 1
#define RESULT_SIZE 256

__thread_local volatile unsigned long get_reply, put_reply, buffer_reply;
__thread_local volatile unsigned long start,end;
__thread_local int my_id; 
__thread_local unsigned long nn_time[5];


static inline unsigned long rpcc()
{
	unsigned long time;
	asm("rtc %0": "=r" (time) : );
	return time;
}
struct Photon
{
	float pos[3];
	unsigned int right;
	unsigned char data[8];
	unsigned short dep;
	char flag;
};
struct Intersection
{
	float p[3];
	float n[3];
};
struct SearchResult
{
	int index;
	float distSquared;
};
struct Morton_node
{
	int code;
	int range[2]; //st ed
};
typedef struct
{
	struct Morton_node *M;
	int *M_S;
	int *T;
	int *S_D;
	int *S_N;
} Info_sort;

extern float m_cosTheta[256];
extern float m_sinTheta[256];
extern float m_cosPhi[256];
extern float m_sinPhi[256];
extern float m_expTable[256];

extern float m_scale;
extern float radius;
extern int task_num[5];
extern float *cut_aabb;
extern int *id_list;

extern float *col_list_task;
extern struct Photon  *sub_tree;
extern struct Intersection *task_its_list;
extern unsigned long long size;

struct Photon* sub_tree_s;
struct Intersection *its_list_s;

__thread_local struct Intersection its_slave[8];
__thread_local float x_exp[11] = {0.03125, 0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16};


__thread_local float s_cosTheta[256];
__thread_local float s_sinTheta[256];
__thread_local float s_cosPhi[256];
__thread_local float s_sinPhi[256];
__thread_local float col[4];

__thread_local float cut_aabb_s[6*300];
__thread_local int id_list_s[300];
__thread_local cut_num_s;
enum {
	ELeafFlag = 0x10,
	EAxisMask = 0x0F
};
inline void getDirection(unsigned char *data, float *res) {
	res[0] = s_cosPhi[data[5]] * s_sinTheta[data[4]];
	res[1] = s_sinPhi[data[5]] * s_sinTheta[data[4]];
	res[2] = s_cosTheta[data[4]];
}
inline void getNormal(unsigned char *data, float *res) {
	res[0] = s_cosPhi[data[7]] * s_sinTheta[data[6]];
	res[1] = s_sinPhi[data[7]] * s_sinTheta[data[6]];
	res[2] = s_cosTheta[data[6]];
}
inline float dot_vv(float *v1, float *v2){
	return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}
inline void dot_vf(float *v1, float t, float *res){
	res[0] = v1[0] * t;
	res[1] = v1[1] * t;
	res[2] = v1[2] * t;
}

inline void sub(float *v1, float *v2, float *res){
	res[0] = v1[0] - v2[0];
	res[1] = v1[1] - v2[1];
	res[2] = v1[2] - v2[2];
}
float my_pow(int y){
	int i = 0;
	float res = 1.0;
	if(y>-6&&y<5){
		return x_exp[y+5];
	}
	if(y>=0){
		for(i;i<y;i++){
			res*=2.0f;
		}
	}
	else{
		float inv_x = 0.5f;
		int y = -y;
		for(i = 0;i<y;i++){
			res*=inv_x;
		}
	}
	
	return res;
} 

void getPower(unsigned char *rgbe, float *result)
{
	if (rgbe[3]) {
		/* Calculate exponent/256 */
		float exp = my_pow((int)rgbe[3] - 136);
		result[0] = rgbe[0] * exp;  //可能有问题 随便找个4096
		result[1] = rgbe[1] * exp;
		result[2] = rgbe[2] * exp;
	}
	else {
		result[0] = result[1] = result[2] = 0.0f;
	}
}
int pos_aabb_intsect(float *pos, float *aabb){   //1包围2 返回1 否则返回0
	int i;
	for (i = 0; i < 3; i++){
		if (pos[i] < aabb[i] || pos[i] > aabb[3 + i])
			return 0;
	}
	return 1;
}
inline int isLeaf(char flags)  { return flags & (unsigned char)ELeafFlag; }
inline unsigned short getAxis(char flags)  { return flags & (unsigned char)EAxisMask; }

void nnSearch(struct Intersection *its_list, float *_sqrSearchRadius,
	size_t k, struct SearchResult *results,struct Photon *photon_tree, int *resultCount) {
	
	int stack[1000];
	char mask_stack[1000];
	int index = 0, stackPos = 1;
	float sqrSearchRadius[8];// = _sqrSearchRadius[i];
	float radius_Max[8];
	int radius_Max_index[8];
	stack[0] = 0;
	struct Photon node; 
	mask_stack[0] = 0x15;
	char mask,child_mask[2];
	char N = 0x15;
	int i,t;
	char T[] = { 0x1, 0x2, 0x4, 0x8 };
	for (i = 0; i < 8; i++){ 
		sqrSearchRadius[i] = _sqrSearchRadius[i]; 
		radius_Max[i] = 0;
		radius_Max_index[i] = 0;
	}
	
	struct Photon buffer[2];
	char mask_buffer[2];
	int index_buffer[2];
	
	int loop = 0;
	int is_buffer = 0;
	for ( t = 0; t<cut_num_s; t++){
		if(pos_aabb_intsect(its_list[0].p, &cut_aabb_s[6*(cut_num_s - t-1)])==0)
			continue;
		
		stack[0] = id_list_s[t];//0;//id_list_s[t];//id_list
		stackPos = 1;
		mask_stack[0] = 0x15;
		
		while (stackPos > 0|| is_buffer!=0) {
			if (is_buffer == 0){
				index = stack[--stackPos];
				
				get_reply = 0;
				athread_get(PE_MODE,&(photon_tree[index]),&(node),sizeof(struct Photon),&get_reply,0,0,0);
				while(get_reply!=1);
			
				mask = mask_stack[stackPos];

			}
			else{
				while(buffer_reply!=1);
				index = index_buffer[loop % 2];
				node = buffer[loop%2];			
				mask = mask_buffer[loop % 2];
				loop++;
				is_buffer = 0;
			}

			if (stackPos > 0){   //write buffer
				index_buffer[loop % 2] = stack[--stackPos];
				
				buffer_reply = 0;
				athread_get(PE_MODE,&(photon_tree[index_buffer[loop % 2]]),&(buffer[loop % 2]),sizeof(struct Photon),&buffer_reply,0,0,0);
				
				mask_buffer[loop % 2] = mask_stack[stackPos];
				is_buffer = 1;
			}	
			
			child_mask[0] = 0x0; child_mask[1] = 0x0;

			for (i = 0; i < PACKET_SIZE; i++){
				if (!mask&N)
					continue;
				struct Intersection *its = &its_list[i];
				if (!isLeaf(node.flag)) {
					float distToPlane = its->p[getAxis(node.flag)] - node.pos[getAxis(node.flag)];
					int searchBoth = distToPlane*distToPlane <= sqrSearchRadius[i];
					if (distToPlane > 0) {
						if(node.right != 0){
							if(searchBoth)
								child_mask[0] |= T[i];
							child_mask[1] |= T[i];
						}
						else if(searchBoth){
							child_mask[0] |= T[i];
						}
								
					}
					else{
						if (searchBoth && node.right != 0)
							child_mask[1] |= T[i];
						child_mask[0] |= T[i];
					}
				}
				
				/* Check if the current point is within the query's search radius */
				float t[3];
				sub(node.pos, its->p, t);
				float pointDistSquared = dot_vv(t, t);

				if (pointDistSquared < sqrSearchRadius[i]) {
					/* Switch to a max-heap when the available search
					result space is exhausted */
					if (resultCount[i] < k) {
						/* There is still room, just add the point to
						the search result list */
						results[i*k + resultCount[i]].distSquared = pointDistSquared;
						if (pointDistSquared > radius_Max[i]){
							radius_Max[i] = pointDistSquared;
							radius_Max_index[i] = resultCount[i];
						}

						results[i*k + resultCount[i]++].index = index;
					}
					else {
						/* Add the new point, remove the one that is farthest away */
						int min_index = 0;
						if (pointDistSquared < radius_Max[i]){
							results[i*k + radius_Max_index[i]].distSquared = pointDistSquared;
							results[i*k + radius_Max_index[i]].index = index;
							radius_Max_index[i] = 0;
							radius_Max[i] = 0;
							int j;
							for (j = 0; j < k; j++){
								if (results[i*k + j].distSquared>radius_Max[i]){
									radius_Max[i] = results[i*k + j].distSquared;
									radius_Max_index[i] = j;
								}
							}
						}

						/* Reduce the search radius accordingly  应该是最近的距离 */
						sqrSearchRadius[i] = radius_Max[i];
					}
				}
			}	
			///处理mask stack  以及stack
			if (child_mask[0] > 0){
				mask_stack[stackPos] = child_mask[0];
				stack[stackPos++] = index + 1;
			}
				
			if (child_mask[1] > 0){
				mask_stack[stackPos] = child_mask[1];
				stack[stackPos++] = node.right;		
			}
				
		}
	
	}
	for (i = 0; i < PACKET_SIZE; i++){
		_sqrSearchRadius[i] = sqrSearchRadius[i];
	}
}
void estimateIrradiance(
	float m_scale, struct Intersection *its_list, 
	float searchRadius, int maxDepth,
	int maxPhotons,struct Photon *photon_tree,float *res,int type) {

	struct SearchResult results[RESULT_SIZE];
	float squaredRadius[8];
	int i,j;
	for (i = 0; i < PACKET_SIZE; i++){
		squaredRadius[i] = searchRadius*searchRadius;
	}
	
	int resultCount[8];
	for(i = 0;i<PACKET_SIZE;i++)
		resultCount[i] = 0;
	nnSearch(its_list, squaredRadius, maxPhotons, results, photon_tree, resultCount);
	

	for(i = 0;i<PACKET_SIZE;i++){
		float invSquaredRadius = 1.0f / squaredRadius[i];
		struct Intersection *its = &its_list[i];
		float *n = its->n;
		float result[] = { 0.0f, 0.0f, 0.0f };
		struct Photon photon;
		for (j = 0; j<resultCount[i]; j++) {
			struct SearchResult *searchResult =  &results[j + i*maxPhotons];
			get_reply = 0;
			athread_get(PE_MODE,&(photon_tree[searchResult->index]),&(photon),sizeof(struct Photon),&get_reply,0,0,0);
			while(get_reply!=1);

			if (photon.dep > maxDepth)
				continue;

			float wi[3], photonNormal[3];
			getDirection(photon.data, wi);
			getNormal(photon.data, photonNormal);
			
			wi[0] = -wi[0]; wi[1] = -wi[1]; wi[2] = -wi[2];
			float wiDotGeoN = fabs(dot_vv(photonNormal, wi)),
				wiDotShN = dot_vv(n, wi);

			if (dot_vv(photonNormal, n) > 1e-1f && wiDotGeoN > 1e-2f) {
				float power[3];
				getPower(photon.data, power);
				float temp = fabs(wiDotShN / wiDotGeoN);
				dot_vf(power, temp, power);

				float sqrTerm = 1.0f - searchResult->distSquared*invSquaredRadius;
				int k;
				for (k = 0; k<3; k++){
					result[k] += power[k] * (sqrTerm*sqrTerm);
				}
			}
		}
		if(type==2){
//			dot_vf(result, m_scale * 3 * INV_PI * invSquaredRadius, &res[0]);
			res[i*4] = result[0]; res[i*4+1] = result[1]; res[i*4+2] = result[2];
			res[i*4+3] = squaredRadius[i];
		}
		else{
			dot_vf(result, m_scale * 3 * INV_PI * invSquaredRadius, &res[i*3]);
		}
		
	}
}
void NNsearch(int *type)
{	
	int i,j;
	get_reply = 0;
	int r_size = RESULT_SIZE;
	sub_tree_s = ( struct Photon*)sub_tree; 
	its_list_s = ( struct Intersection*)task_its_list; 
	
	float *col_list_task_s = (float*)col_list_task;
	int task_num_s;

	task_num_s = task_num[0];
	cut_num_s = task_num[2];
	float m_scale_s = m_scale;//9.64348055e-006;
	float radius_s = radius;

	my_id = athread_get_id(-1);
	int PS = PACKET_SIZE;
	
	get_reply = 0;
	athread_get(PE_MODE, &(m_cosTheta[0]), &(s_cosTheta[0]), 256*4, &get_reply, 0, 0, 0);
	athread_get(PE_MODE, &(m_sinTheta[0]), &(s_sinTheta[0]), 256*4, &get_reply, 0, 0, 0);
	athread_get(PE_MODE, &(m_cosPhi[0]), &(s_cosPhi[0]), 256*4, &get_reply, 0, 0, 0);
	athread_get(PE_MODE, &(m_sinPhi[0]), &(s_sinPhi[0]), 256*4, &get_reply, 0, 0, 0);
	athread_get(PE_MODE, &(cut_aabb[0]), &(cut_aabb_s[0]), cut_num_s*6*4, &get_reply, 0, 0, 0);
	athread_get(PE_MODE, &(id_list[0]), &(id_list_s[0]), cut_num_s*4, &get_reply, 0, 0, 0);
	while(get_reply!=6);
	

	for (i = my_id*PS;i < task_num_s; i = i + 64*PS) {
		get_reply = 0;
		athread_get(PE_MODE, &(its_list_s[i]), &(its_slave[0]), 24*PS, &get_reply, 0, 0, 0);
		while(get_reply!=1);
		
//		estimateIrradiance(m_scale_s, &its_slave[0], radius_s, 20, r_size, sub_tree_s, col, *type);
		
		put_reply = 0;
		if(*type==2)
			athread_put(PE_MODE, &(col[0]), &(col_list_task_s[i*4]), 16*PS, &put_reply, 0, 0);
		else
			athread_put(PE_MODE, &(col[0]), &(col_list_task_s[i*3]), 12*PS, &put_reply, 0, 0);
		//可以等到最后一起put
		while(put_reply!=1); 
		
	}	
}
int partition(struct Morton_node *list, int top, int bottom)
{
	int x = list[top].code;
	int i = top - 1;
	int j = bottom + 1;
	int temp;
	do
	{
		do
		{
			j--;
		} while (x < list[j].code);
		do
		{
			i++;
		} while (x >list[i].code);

		if (i < j)
		{
			struct Morton_node temp = list[i];
			list[i] = list[j];
			list[j] = temp;

		}
	} while (i < j);
	return j;           // returns middle subscript  
}

void quickSort(struct Morton_node *list, int top, int bottom)
{
	// top = subscript of beginning of array
	// bottom = subscript of end of array

	int middle;
	if (top < bottom)
	{
		middle = partition(list, top, bottom);
		quickSort(list, top, middle);   // sort first section
		quickSort(list, middle + 1, bottom);    // sort second section
	}
	return;
}
void PSRS_sort_part1(struct Morton_node *morton_list){
	my_id = athread_get_id(-1);
	int sort_task_num = size / 64;
	int sample_num = sort_task_num / 64;
	int i,j;

	quickSort(&morton_list[my_id*sort_task_num], 0, sort_task_num - 1);	
	
	struct Morton_node *main_element = &morton_list[size+64*my_id];
	for (j = 0; j < 64; j++){
		main_element[j].code = morton_list[my_id*sort_task_num + j*sample_num].code;  
	}
}
void PSRS_sort_part2(Info_sort *_ptr){
	my_id = athread_get_id(-1);
	int sort_task_num = size / 64;
	int sample_num = sort_task_num / 64;
	
	Info_sort *info = _ptr;	
	struct Morton_node *morton_list = &info->M[0];
	int *main_sample = info->M_S;
	int *tag = &info->T[my_id*65];
	int date_range = sort_task_num*2;
	int *sorted_data = &info->S_D[my_id * date_range];	
	int *thread_tag = &info->T[my_id];
	int current = 0;
	int j = 0;
	tag[0] = my_id*sort_task_num;
	for (j = 1; j < 65; j++){
		int current_main_tag = main_sample[j];    //小于这个tag的
		while (morton_list[my_id*sort_task_num+current].code < current_main_tag && current < sort_task_num){
			current++;
		}
		tag[j] = current + my_id*sort_task_num;
	}

}
















