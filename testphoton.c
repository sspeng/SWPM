#include "stdio.h"
#include "math.h"
#include <athread.h>
#include<string.h>
#include "tree.h"
#include "cluster.h"
#include "mpi.h"
#define M_PI 3.14159265358979323846
#define INV_PI 0.31830988618379067154
#define inline __inline
////only athread 
#define MASTER_WORK 0   //
#define OVERLAP_ALL_PHOTON 1     //overlap all photon or only viewpoint direction
#define MSG_TAG 10
#define UNIT 15000000
#define INT_MAX 2147483647
#define FLOAT_MAX 3.402823e+038
#define FLOAT_MIN 1.175494e-038
#define BVH_TREE_LV 14
#define MAX_ITS_LEAF_NUM 256
//#define USE_SLAVE 4
float m_cosTheta[256];
float m_sinTheta[256];
float m_cosPhi[256];
float m_sinPhi[256];
float m_expTable[256];

extern SLAVE_FUN(NNsearch)();
extern SLAVE_FUN(NNsearch_tree)();
static inline unsigned long rpcc()
{
	unsigned long time;
	asm("rtc %0": "=r" (time) : );
	return time;
}

float *cut_aabb;
int *id_list;
struct Intersection  *its_list;
struct Intersection *task_its_list;
int *task_its_index;
struct Photon *x_tree;
struct Photon *sub_tree;
float *col_list;//[3*512*512];
float *col_list_task; 
int task_num[5] = {0};    //0 task num  1 sub_tree size 
float m_scale, radius = 0;
int its_num = 1024*1024;   //task number
unsigned long long size, depth;
unsigned long st,ed;
unsigned long all_st, all_ed, unit;
unsigned long slave_work_time = 0;
unsigned long time_statistics[3];
unsigned long data_statistics[3];
int mpi_id, numprocs, slave_num;
MPI_Request request[10];
MPI_Group G_psrs;   //for psrs
MPI_Comm comm_psrs;
MPI_Status Istatus[10] ;
struct SearchResult
{
	int index;
	float distSquared;
}results[1000];

void createPrecompTables() {
	int i;
	for (i = 0; i<256; i++) {
		double angle = (float)i * ((float)M_PI / 256.0f);
		m_cosPhi[i] = cos(2.0f * angle);
		m_sinPhi[i] = sin(2.0f * angle);
		m_cosTheta[i] = cos(angle);
		m_sinTheta[i] = sin(angle);
		m_expTable[i] = pow(2.0, i - (128 + 8));
/*	*/
	}
	m_expTable[0] = 0;

}
inline void getDirection(unsigned char *data, float *res) {
	res[0] = m_cosPhi[data[5]] * m_sinTheta[data[4]];
	res[1] = m_sinPhi[data[5]] * m_sinTheta[data[4]];
	res[2] = m_cosTheta[data[4]];
}
inline void getNormal(unsigned char *data, float *res) {
	res[0] = m_cosPhi[data[7]] * m_sinTheta[data[6]];
	res[1] = m_sinPhi[data[7]] * m_sinTheta[data[6]];
	res[2] = m_cosTheta[data[6]];
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
inline void add(float *v1, float *v2, float *res){
	res[0] = v1[0] + v2[0];
	res[1] = v1[1] + v2[1];
	res[2] = v1[2] + v2[2];
}
void getPower(unsigned char *rgbe, float *result)
{
	if (rgbe[3]) {
		/* Calculate exponent/256 */
		float exp = pow(2.0, (int)rgbe[3] - (128 + 8));
		result[0] = rgbe[0] * exp;  //可能有问题 随便找个4096
		result[1] = rgbe[1] * exp;
		result[2] = rgbe[2] * exp;
	}
	else {
		result[0] = result[1] = result[2] = 0.0f;
	}
}
void build_struct_photon(MPI_Datatype *mpi_photon){
	MPI_Datatype type[5] = { MPI_FLOAT, MPI_UNSIGNED, MPI_UNSIGNED_CHAR, MPI_UNSIGNED_SHORT,MPI_CHAR};
	int blocklen[5] = { 3, 1, 8, 1, 1 };
    MPI_Aint address[5];
	MPI_Aint disp[5];
	
	MPI_Get_address(&x_tree[0].pos[0],&address[0]);
	MPI_Get_address(&x_tree[0].right,&address[1]);
	MPI_Get_address(&x_tree[0].data[0],&address[2]);
	MPI_Get_address(&x_tree[0].dep,&address[3]);
	MPI_Get_address(&x_tree[0].flag,&address[4]);
	int i = 0;
	for(i;i<5;i++){
		disp[i] = address[i] - address[0];
	}
	MPI_Type_struct(5,blocklen,disp, type, mpi_photon);
	MPI_Type_commit(mpi_photon);
}
void build_struct_bvh(MPI_Datatype *mpi_bvh_tree,struct BVH_tree *node){
	MPI_Datatype type[2] = { MPI_FLOAT, MPI_INT};
	int blocklen[2] = { 6, 5 };
    MPI_Aint address[2];
	MPI_Aint disp[2];
	
	MPI_Get_address(&node[0].box[0],&address[0]);
	MPI_Get_address(&node[0].children[0],&address[1]);

	int i = 0;
	for(i;i<2;i++){
		disp[i] = address[i] - address[0];
	}
	MPI_Type_struct(2,blocklen,disp, type, mpi_bvh_tree);
	MPI_Type_commit(mpi_bvh_tree);
}
int in_sub_tree(float aabb_box[6], float pos[3]){
	int k;
	for (k = 0; k < 3; k++){
		if (pos[k]<aabb_box[k] || pos[k] > aabb_box[k+3])
			return 0;
	}
	return 1;
}
void sub_bound_box(float min[3], float max[3], float *list, int *num, int iterator){
	if (iterator == -1) return;
	if (iterator == 0){
		list[(*num) * 6] = min[0];
		list[(*num) * 6 + 1] = min[1];
		list[(*num) * 6 + 2] = min[2];
		list[(*num) * 6 + 3] = max[0];
		list[(*num) * 6 + 4] = max[1];
		list[(*num) * 6 + 5] = max[2];
		(*num)++;
		return;
	}
	int axis = getLargestAxis_old(max, min);
	float midpoint[3];
	midpoint[axis] = 0.5f * (max[axis] + min[axis]);
	int k = 0;
	for (k = 0; k < 3; k++){
		if (k != axis)
			midpoint[k] = max[k];
	}
	sub_bound_box(min, midpoint,  list, num, iterator - 1);
	for (k = 0; k < 3; k++){
		if (k != axis)
			midpoint[k] = min[k];
	}
	sub_bound_box(midpoint, max, list, num, iterator - 1);
}
void read_photon_data(char *path, float *aabb_box){
	char photon_path[100];
	strcpy(photon_path,path);
	strcat(photon_path,"photonmap");  
	FILE *photon_input = fopen(photon_path, "rb");
		
	fread(&m_scale, sizeof(float), 1, photon_input);
	fread(&size, sizeof(unsigned long long), 1, photon_input);
	fread(&depth, sizeof(unsigned long long), 1, photon_input);
	fread(aabb_box, sizeof(float), 6, photon_input);
		
	x_tree=(struct Photon*)malloc((size+100)*sizeof(struct Photon));
	int i;
	for (i = 0; i < size; i++){
		fread(&x_tree[i].pos[0], sizeof(float), 3, photon_input);
		fread(&x_tree[i].right, sizeof(unsigned int), 1, photon_input);
		fread(&x_tree[i].data[0], sizeof(unsigned char), 8, photon_input);
		fread(&x_tree[i].dep, sizeof(unsigned short), 1, photon_input);
		fread(&x_tree[i].flag, sizeof(char), 1, photon_input);
	}
	fclose(photon_input);
}
void read_its_data(char *path){
	char its_path[100];
	strcpy(its_path,path);
	strcat(its_path,"its");  
	FILE *its_input = fopen(its_path, "rb");
	fread(&radius, sizeof(float), 1, its_input);
	fread(&its_num, sizeof(int), 1, its_input);
	
	its_list=(struct Intersection*)malloc((its_num) * sizeof(struct Intersection));
	int num = 0;
	while(!feof(its_input)){
		fread(&its_list[num].p[0], sizeof(float), 3, its_input);
		fread(&its_list[num].n[0], sizeof(float), 3, its_input);
		num++;
	}
//	printf("%f %f %f \n",its_list[0].p[0],its_list[0].p[1],its_list[0].p[2]);
//	printf("old its %d new its %d list_size %d\n",its_num, num, its_num*2);
	its_num = num - 1;
	fclose(its_input);
}
void free_all(){
	free(x_tree);
	free(its_list);
	free(task_its_list);
	free(task_its_index);
	free(col_list);
	free(col_list_task);
}

void normalize_color(){
	int i = 0, k = 0;
	float max_col[3] = {0}, min_col[3] = {0};
	for(i = 0;i<its_num;i++){
//		printf("%d max_col %f %f %f, min_col %f %f %f \n", i,max_col[0],max_col[1],max_col[2],min_col[0],min_col[1],min_col[2]);
		for(k = 0;k<3;k++){
			if(col_list[i*3+k] < min_col[k]||i == 0)
				min_col[k] = col_list[i*3+k];
			if(col_list[i*3+k] > max_col[k]||i == 0)
				max_col[k] = col_list[i*3+k];
		}
	}
	printf("max_col %f %f %f, min_col %f %f %f \n", max_col[0],max_col[1],max_col[2],min_col[0],min_col[1],min_col[2]);
	float range = 0;
	for(k = 0;k<3;k++){
		float temp = max_col[k] - min_col[k];
		if(temp>range)
			range = temp;
	}
	float insrange = 1/range;
	printf("max range %f \n",range);
	for(i = 0;i<its_num;i++){
		for(k = 0; k < 3; k++){
			col_list[i*3 + k] *= insrange;
		}
	}
}
void write_color(char *path){
	char reslut_path[100];
	strcpy(reslut_path, path);
	strcat(reslut_path,"result");  
	FILE *resoutput = fopen(reslut_path,"wb");
	int i;
	for(i = 0;i<its_num;i++){
		fwrite(&col_list[i*3], sizeof(float), 3, resoutput);
	}
	fclose(resoutput);
}
void print_static(double statistics[5][256], double time_static[10], 
		long long all_photon, int numprocs, int mul, char *path, int type){
			
	int j,i;
	double average[5], variance[5];
	for(i = 0;i<5;i++){average[i] = 0; variance[i] = 0;}
	for(i = 0;i< numprocs-1;i++){
		for(j = 0;j<5;j++){
			average[j] += statistics[j][i];
		}
	}
	for(j = 0;j<5;j++){
		average[j] /= numprocs-1;
	}
	printf("slave averge tree construction time %f\n",average[0]);
	for(i = 0;i< numprocs-1;i++){
		for(j = 0;j<5;j++){
			variance[j] += (average[j]-statistics[j][i])*(average[j]-statistics[j][i]);
		}
	}
	for(j = 0;j<5;j++){
		variance[j] /= numprocs-1;
	}
	
	char reslut_path[100], str[25];
	strcpy(reslut_path, path);
	strcat(reslut_path,"timestatic");
	
	FILE *access = fopen(reslut_path,"r");
	FILE *resoutput = fopen(reslut_path,"a");
	if(!access){				
		fprintf(resoutput, "preprocess | dis prepare |  dis   |  wait  | master time || slave time | bulid tree | estimate | photon size | its size\n");
	}
	else{
		fclose(access);
	}
	if(type==0)
		fprintf(resoutput, "overlap slave: %d  task granularity: %d  allphoton %ld\n", numprocs, mul, all_photon);
	else
		fprintf(resoutput, "tree  slave: %d  task granularity: %d  allphoton %ld\n", numprocs, mul, all_photon);
	fprintf(resoutput, "%5.2f          %5.2f       %5.2f     %5.2f    %5.2f        %5.2f      %5.2f     %5.2f    %5.2e     %5.2e \n",
					time_static[0], time_static[1], time_static[2], time_static[3], time_static[4], average[0],average[1],average[2],average[3],average[4]);
	fprintf(resoutput, "                                                            %5.2f      %5.2f     %5.2f    %5.2e     %5.2e \n", 
					variance[0],variance[1],variance[2],variance[3],variance[4]);
	
	fclose(resoutput);
}

float large_radius(float *aabb){
	int i = 0;
	float res = 0;
	for(i;i<3;i++){
		float temp = aabb[3+i]-aabb[i];
		if(temp>res){
			res = temp;
		}
	}
	return res;
}
float min_radius(float *aabb){
	int i = 0;
	float res = 0;
	for(i;i<3;i++){
		float temp = aabb[3+i]-aabb[i];
		if(temp<res){
			res = temp;
		}
	}
	return res;
}
void div_its_aabb(int ab_head, int ab_rear, struct Morton_node *morton_list,
			int *its_num_queue, float *aabb_queue, int *range_queue){		
	//包围盒分割
	float *head_aabb = &aabb_queue[ab_head * 6];
	float *rear_aabb = &aabb_queue[ab_rear * 6];
	
	int i = 0, div_id;
	div_id = getLargestAxis_old(&head_aabb[3], head_aabb);
	for(i = 0; i < 6; i++)
	{
		rear_aabb[i] = head_aabb[i];
		printf("%f",head_aabb[i]);
		
	}
	float mid = (head_aabb[div_id + 3] + head_aabb[div_id]) / 2;
	head_aabb[div_id + 3] = mid;
	rear_aabb[div_id] = mid;

	
	for(i = 0; i<6; i++)
	{
		printf("%f ",head_aabb[i]);
	}
//	printf("head \n");
	for(i = 0; i<6; i++)
	{
		printf("%f ",rear_aabb[i]);
	}
//	printf("rear \n");
	
//	printf("old head num %d \n",its_num_queue[ab_head]);
	
	
	its_num_queue[ab_rear] = its_num_queue[ab_head];
	
	
	
/*	int st = range_queue[ab_head*2], ed = range_queue[ab_head*2+1];
//	printf("range %d  %d \n",st,ed);	
	int temp_num[2] = {0};
	for (i = st; i < ed; i++){
		struct Intersection *its_node = &its_list[morton_list[i].its_id];
		int j = 0;
		if(pos_aabb_intsect(its_node->p, head_aabb)){
			temp_num[0]++;
		}else{
			//(pos_aabb_intsect(its_list[i].p, rear_aabb)){
			temp_num[1]++;
		}
	}
	printf("temp_num[0] %d  temp_num[1]%d  \n",temp_num[0],temp_num[1]);
	if(temp_num[0]==0&&temp_num[1]==0){
		its_num_queue[ab_rear] = 0;
		head_aabb[div_id + 3] = rear_aabb[div_id+3];
		printf("all 0, %d \n", its_num_queue[ab_head]);
	}
	else if(temp_num[0]<5||temp_num[1]<5){
		printf("rear set 0\n");
		its_num_queue[ab_rear] = 0;
		its_num_queue[ab_head] = temp_num[1]+temp_num[0];
	}else{
		its_num_queue[ab_head] = temp_num[0];
		its_num_queue[ab_rear] = temp_num[1];
	}*/
	
	printf("%d head num %d, %drear num %d \n",ab_head,its_num_queue[ab_head],ab_rear, its_num_queue[ab_rear]);
}

void master_pre_estimation(struct Morton_node *morton_list, float *aabb_queue,int *its_num_queue,
					int task_queue_num,	struct BVH_tree *photon_bvh_tree, char *path){
	
	MPI_Status status ; 
	MPI_Status Istatus[10] ; 
	printf("read \n");
	
	unsigned long wait_st,wait_ed,wait_total = 0,data_dis = 0, all_st, all_ed;
	long long all_photon = 0;
	double statistics[5][256], time_static[10];
	int max_its = 0, max_photon = 0, min_its = its_num, min_photon = size;
	
	int k, i,j;

//	printf("read over radius %f\n", radius);
	all_st = rpcc();
	st=rpcc();
	int *power = (int*)malloc(sizeof(int) * 3) ; memset(power,0,sizeof(int)*3) ;  	
	int *its_id_temp = (int*)malloc(its_num*sizeof(int));
	for(i = 0;i<its_num;i++) 
		its_id_temp[i] = i;	
	
	double two = 2;
	
	col_list = (float*)malloc(its_num*3*sizeof(float));	
	memset( col_list, 0, its_num*3*sizeof(float));
	
	float *col_buffer = (float*)malloc(its_num*2*3*sizeof(float));	
	
	int task_num_buffer[2];                                             //double buffer
	int buffer_size = its_num/2;
		
	ed=rpcc();
	time_static[0] = (ed-st)/unit;time_static[0]/=100;
	
	int rece_num = 0, buffer_st = 0, its_task_num = 0;
	st=rpcc();
	int cut_num = 0, subtree_num = 0;
	int top_tree_size = pow(2.0, BVH_TREE_LV) - 1;
	int *cut = (int*)malloc(top_tree_size*6*sizeof(int));  //  cut 直接保存 要读取的文件名 以及数据的索引就行
	cut_aabb = (float*)malloc(top_tree_size*12*sizeof(float));
	int max_subtree = 150000000;//((size*4)/mul)/(numprocs-1);//size;//
	int ab_rear = task_queue_num, ab_head = 0, range_id = 0;
	int max_its_num = 4*its_num/(numprocs-1);
	
	
	float all_aabb[6];
	all_aabb[0] = -1000;all_aabb[1] = -1000;all_aabb[2] = -1000;
	all_aabb[3] = 1000;all_aabb[4] = 1000;all_aabb[5] = 1000;
	
	int *task_range_queue = (int*)malloc(task_queue_num*4*sizeof(int));
	
	//  task range;
	task_range_queue[0] = 0;
	task_range_queue[1] = its_num_queue[0];
	for(i = 1;i<task_queue_num;i++){
		task_range_queue[i*2] = task_range_queue[2*i-1];
		task_range_queue[i*2+1] = task_range_queue[i*2]+its_num_queue[i];
//		printf("task range %d  %d \n",task_range_queue[i*2],task_range_queue[i*2+1]);
	}

	while(ab_head <= ab_rear){
//		printf(" head %d rear %d \n", ab_head, ab_rear);
		if(its_num_queue[ab_head]==0){
			ab_head++;
			continue;
		}
/*		if(its_num_queue[ab_head]>max_its_num){
			
			div_its_aabb(ab_head, ab_rear,morton_list, its_num_queue, aabb_queue, task_range_queue);	
			if(task_range_queue[ab_head*2]>0)   
				task_range_queue[ab_head*2] = -task_range_queue[ab_head*2];   //分割过的是负值					
			if(its_num_queue[ab_rear]!=0){
				task_range_queue[ab_rear*2] = -task_range_queue[ab_rear*2];
				task_range_queue[ab_rear*2+1] = -task_range_queue[ab_rear*2+1];			
				ab_rear++;
			}	
			printf("devid 1 \n");
		}	*/
		float *node_aabb = &aabb_queue[ab_head*6];
//		float *node_aabb = &all_aabb[0];
		
		int *node_range = &task_range_queue[ab_head*2];
				
		float current_radius = radius;
		do{
			cut_num = 0;
			subtree_num = 0;										
			float big_node_aabb[6] = {
				node_aabb[0] - current_radius, node_aabb[1] - current_radius, node_aabb[2] - current_radius ,
				node_aabb[3] + current_radius, node_aabb[4] + current_radius, node_aabb[5] + current_radius  };
			
			cut_num = traverse_BVH_tree(photon_bvh_tree, big_node_aabb, cut, cut_aabb, radius, top_tree_size);
							//检查各个文件是否正确   如果正确则是遍历   建树的问题 
			for (i = 0; i < cut_num; i++){					
				subtree_num += (cut[i*3+2] - cut[i*3+1]);
			}		

/////////////此处很可能有错误			
			if(subtree_num>max_subtree||cut_num>300){   //
		//		float large_axis = large_radius(node_aabb);
				float min_axis = large_radius(node_aabb);
				if(min_axis>2*radius){					
					if(task_range_queue[ab_head*2]<0){
						printf("negetive\n");
						task_range_queue[ab_head*2] = -task_range_queue[ab_head*2];	
					}
					// task_range_queue[ab_head] is real value
					div_its_aabb(ab_head, ab_rear,morton_list, its_num_queue, aabb_queue, task_range_queue);
					   
					task_range_queue[ab_head*2] = -task_range_queue[ab_head*2];   //分割过的是负值						
					task_range_queue[ab_rear*2] = -task_range_queue[ab_rear*2];
					task_range_queue[ab_rear*2+1] = task_range_queue[ab_rear*2+1];
					printf("divide2,cut_num %d subtree_num %d\n",cut_num,subtree_num);
					ab_rear++;			
				}
				else{
					current_radius/=2;
					printf("old %f cur %f subtree_num %d\n", radius, current_radius, subtree_num);
				}
				continue;
			}
//			printf("subtree_num %d \n",subtree_num);
		}while(subtree_num>max_subtree&&current_radius>0.00001);//
		ab_head++;
		
		if(current_radius<0.00001)
			printf("radius too small\n");
		all_photon += subtree_num;
		if(cut_num>300||cut_num==0){
			printf("cut num error %d",cut_num);
			continue;
		}

		task_num[0] = node_range[0];
		task_num[1] = subtree_num;
		task_num[2] = cut_num;
		task_num[3] = node_range[1];
//		printf("its_task_num task_num[0] %d\n",task_num[0]);
			
		max_photon = subtree_num>max_photon?subtree_num:max_photon;
		min_photon = subtree_num<min_photon?subtree_num:min_photon;
						
		wait_st = rpcc();
			//receive i am idle
		MPI_Recv(power,3,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status) ;
									
		MPI_Send(&task_num[0], 4, MPI_INT, power[1], MSG_TAG, MPI_COMM_WORLD);
		wait_ed = rpcc();
		wait_total += (wait_ed - wait_st);
		wait_st = rpcc();
		MPI_Send(&cut[0], cut_num*3, MPI_INT, power[1], MSG_TAG, MPI_COMM_WORLD);		
//		printf("master its %d %f %d cut  %d  %d\n", its_task_num, current_radius, power[1], cut_num, subtree_num);
		MPI_Send(&cut_aabb[0], cut_num*6, MPI_FLOAT, power[1], MSG_TAG, MPI_COMM_WORLD);
		MPI_Send(&node_aabb[0], 6, MPI_FLOAT, power[1], MSG_TAG, MPI_COMM_WORLD);
//		printf("master send %d %d\n", power[1], its_num_queue[ab_head-1]);
		wait_ed = rpcc();
		data_dis += (wait_ed - wait_st);

	}
	ed=rpcc();
	time_static[1] = (ed-st-wait_total-data_dis)/unit; time_static[1] /= 100;
	
	printf("Data organization %f\n", time_static[1]);
	
	st=rpcc();
	int terminate = 0;
	task_num[0] = -1;
	
	printf("all send is over\n");
	printf("All wait1 %ld\n  ", wait_total/unit);
	int *task_its_index_buffer = (int*)malloc(its_num*sizeof(int));
	unsigned long slave_read_time = 0, slave_athread_time = 0,slave_wait_time = 0;
	while(terminate<slave_num){
		wait_st = rpcc();
		MPI_Recv(power,3,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status) ;
			
		MPI_Send(&task_num[0], 4, MPI_INT, power[1], MSG_TAG, MPI_COMM_WORLD);	
		
		MPI_Recv(&its_task_num,1,MPI_INT,power[1],MSG_TAG,MPI_COMM_WORLD,&status) ;
		
		max_its = its_task_num>max_its?its_task_num:max_its;
		min_its = its_task_num<min_its?its_task_num:min_its;		
		
		wait_ed = rpcc();
		wait_total += (wait_ed - wait_st);
		wait_st = rpcc();
		MPI_Recv(&task_its_index_buffer[0], its_task_num, MPI_INT, power[1], MSG_TAG, MPI_COMM_WORLD,&status);
		MPI_Recv(&col_buffer[0], its_task_num*3, MPI_FLOAT, power[1], MSG_TAG, MPI_COMM_WORLD,&status);
		for(j = 0;j<its_task_num;j++){
			int new_col_index = task_its_index_buffer[j] * 3;
			col_list[new_col_index] = col_buffer[j*3];
			col_list[new_col_index + 1] = col_buffer[j*3+1];
			col_list[new_col_index + 2] = col_buffer[j*3+2];
		}
		printf("%d    ", its_task_num);
		wait_ed = rpcc();
		data_dis += (wait_ed - wait_st);	
		
		MPI_Recv(&time_statistics[0], 3, MPI_LONG, power[1], MSG_TAG, MPI_COMM_WORLD,&status);
		MPI_Recv(&data_statistics[0], 2, MPI_LONG, power[1], MSG_TAG, MPI_COMM_WORLD,&status);

		slave_read_time += time_statistics[0]; 
		slave_athread_time+= time_statistics[1];
		slave_wait_time += time_statistics[2];
		
		statistics[3][terminate] = data_statistics[0]; 
		statistics[4][terminate] = data_statistics[1];		
		terminate++;
			
	}
	printf("min photon %d its %d, max photon %d its %d \n", min_photon,min_its, max_photon,max_its);
	all_ed = rpcc();
	slave_read_time = slave_read_time/(numprocs-1)/unit;
	slave_athread_time = slave_athread_time/(numprocs-1)/unit;
	slave_wait_time = slave_wait_time/(numprocs-1)/unit;
	
	printf("All wait %ld\n  ", wait_total/unit);
	printf("All Communication %ld\n  ", data_dis/unit);
	printf("slave read time %ld  athread %ld wait time %ld\n",slave_read_time,slave_athread_time,slave_wait_time);
	printf("All used photon %ld\n", all_photon);
	
	free(its_id_temp);
//	free(photon_tree);
	free(power);
}
void cut_read(char *path, struct Photon *tree, int *cut, int cut_num,int *id_list){
	int i,j;
	char photon_path[110];
	strcpy(photon_path, path);
	strcat(photon_path, "new_photon_");
	FILE *photon_input = NULL; 	
	int tree_size = 0;
	id_list[0] = 0;
	cut[cut_num * 3] = -1;
	char file_id[5],read_path[120];
//	printf("cut num %d \n",cut_num);
	for (i = cut_num-1; i >= 0; i--){
//		printf("%d %d %d  ",cut[i * 3],cut[i * 3+1],cut[i * 3+2]);
		if (cut[i * 3] == cut[(i + 1) * 3]){
			fseek(photon_input, (cut[i * 3 + 1] - cut[(i+1) * 3 + 2]) * sizeof(struct Photon), SEEK_CUR);
		}
		else{
			if (i != cut_num-1)
				fclose(photon_input);
			sprintf(file_id, "%d", cut[i*3]);
			strcpy(read_path, photon_path);
			strcat(read_path, file_id);		
			photon_input = fopen(read_path, "rb");
			fseek(photon_input, cut[i * 3 + 1] * sizeof(struct Photon), SEEK_CUR);
		}
		int st = cut[i * 3 + 1];
		int n = cut[i * 3 + 2] - st;
		for (j = 0; j < n; j++){	
			fread(&tree[tree_size], sizeof(struct Photon), 1, photon_input);
			if (!isLeaf(tree[tree_size].flag) && tree[tree_size].right!=0)
				tree[tree_size].right = id_list[cut_num - i-1] + tree[tree_size].right - st;
			tree_size++;
		}
		id_list[cut_num - i] = id_list[cut_num - i - 1] + n;
	}
//	printf("\n");
}

void slave_pre_estimation(struct Morton_node *morton_list, int *queue_range, 
						int local_task_num, struct BVH_tree *photon_bvh_tree, char *path){
	MPI_Status status;  
	int i = 0;
	int *power ; 
	unsigned long read_time = 0, athread_time = 0,wait_time = 0,read_st, read_ed, ath_st, ath_ed;
	int all_photon = 0,all_its = 0;
	power = (int*)malloc(sizeof(int)*3) ;  
	power[0] = 0 ;   power[1] = mpi_id ;  power[2] = 0;
	all_st = rpcc();
	///work  	
	struct Intersection* its_buffer = (struct Intersection*)malloc(its_num*sizeof(struct Intersection));	
	float *col_buffer=(float*)malloc(its_num*3*sizeof(float));
	task_its_index = (int*)malloc(its_num*sizeof(int));
		
	float aabb_its[6] = {0};
	int *its_id_temp = (int*)malloc(its_num*sizeof(int));
	for(i = 0;i<its_num;i++) 
		its_id_temp[i] = i;	
	
	int top_tree_size = pow(2.0, BVH_TREE_LV) - 1;//(numprocs-1)*8-1;
	int *cut = (int*)malloc(top_tree_size*3*sizeof(int));
	float *cut_aabb_buf = (float*)malloc(top_tree_size*12*sizeof(float)); 
	float *its_aabb_buf = (float*)malloc(12*sizeof(float));
	int *id_list_buf = (int*)malloc(top_tree_size*2*sizeof(int));
	int max_tree = size;
	struct Photon* sub_tree_buf;// = (struct Photon*)malloc(max_tree*2*sizeof(struct Photon)); 
	int task_buf[8];
	int loop = 0, buf_st = 0, cut_buf_st;
	int task_its_num;
	
	createPrecompTables();
	int receive_num = 0;
	for(;;){
		st = rpcc();
		power[2] = 1;
		MPI_Send(power, 3, MPI_INT, 0, MSG_TAG,MPI_COMM_WORLD) ; 
		MPI_Recv(&task_buf[buf_st*4], 4, MPI_INT, 0, MSG_TAG, MPI_COMM_WORLD, &status);
		ed = rpcc();
		wait_time+=ed-st;
		if(task_buf[buf_st*4]!=-1){
//			if(local_task_num>0){
//				读取本地数据	
//			}
			read_st = rpcc();
			int subtree_num = task_buf[buf_st*4+1];
			int cut_num = task_buf[buf_st*4+2];
			int range[2];
			range[0] = task_buf[buf_st*4+0];
			range[1] = task_buf[buf_st*4+3];
//			printf("range %d  %d\n",range[0],range[1]);
			
			
			MPI_Recv(&cut[0], cut_num*3, MPI_INT, 0, MSG_TAG, MPI_COMM_WORLD, &status);
			MPI_Recv(&cut_aabb_buf[top_tree_size*6*buf_st], cut_num*6, MPI_FLOAT, 0, MSG_TAG, MPI_COMM_WORLD, &status);
			MPI_Recv(&its_aabb_buf[6*buf_st], 6, MPI_FLOAT, 0, MSG_TAG, MPI_COMM_WORLD, &status);
			
			task_its_num = all_its;
			int test_in_num = 0;
			
/*			for (i = 0; i < its_num; i++){
				if(pos_aabb_intsect(its_list[i].p, &its_aabb_buf[6*buf_st])){
					its_buffer[all_its] = its_list[i];
					task_its_index[all_its++] = its_id_temp[i];	
				}
			}*/
			
			receive_num++;
			if(range[0]<0){
				range[0] = -range[0];
				for (i = range[0]; i < range[1]; i++){
					struct Intersection *its_node = &its_list[morton_list[i].its_id];
//					printf("%d  ",morton_list[i].its_id);
					if(pos_aabb_intsect(its_node->p, &its_aabb_buf[6*buf_st])){
						its_buffer[all_its] = *its_node;
						task_its_index[all_its++] = morton_list[i].its_id;	
					}
				}
			}
			else{
				for (i = range[0]; i < range[1]; i++){
//					printf("%d ",range[0]);
//					printf("%d  ",morton_list[i].its_id);
					struct Intersection *its_node = &its_list[morton_list[i].its_id];					
					its_buffer[all_its] = *its_node;
					task_its_index[all_its++] = morton_list[i].its_id;
				}
			}
			
			
			task_its_num = all_its - task_its_num;
			
			int j;
			
/*			for(j = 0;j<6;j++){
				printf("%f ",its_aabb_buf[6*buf_st + j]);
			} 
			printf("mpi id %d its  %d subtree_num %d \n", mpi_id, task_its_num, subtree_num);
*/			
			if(task_its_num == 0){
				printf(" receive o task\n");
				
				continue;
				
			}
				
			
			
			
			sub_tree_buf = (struct Photon*)malloc(subtree_num*sizeof(struct Photon));
			cut_read(path, sub_tree_buf, cut, cut_num, &id_list_buf[top_tree_size*buf_st]);	
			
			id_list_buf[top_tree_size*buf_st+cut_num] = subtree_num;
//			printf("%d sub tree size %d \n",mpi_id,subtree_num);
			read_ed = rpcc();
			read_time += read_ed - read_st;
			if(loop!=0){
				athread_join();
				ath_ed = rpcc();
				athread_time +=(ath_ed-ath_st);
//				printf("%d athread_join over %d \n",mpi_id,subtree_num);
			}
			task_num[0] = task_its_num; //task_buf[buf_st*4];   
			task_num[2] = task_buf[buf_st*4+2];
			cut_aabb = &cut_aabb_buf[top_tree_size*6*buf_st];
			id_list = &id_list_buf[top_tree_size*buf_st];
			task_its_list = &its_buffer[all_its-task_num[0]];
			col_list_task = &col_buffer[3*(all_its-task_num[0])];
			free(sub_tree);
			sub_tree = sub_tree_buf;
			int type = 1;
			
			ath_st = rpcc();
			athread_spawn(NNsearch,&type); 
			
			buf_st = loop%2;	
		
		}		
		else{	
//			printf("%d athread_join  \n",mpi_id);
//			printf("%d over\n", mpi_id);
			athread_join();	
			ath_ed = rpcc();
			athread_time +=(ath_ed-ath_st);
//			printf("%d final athread_join over %d \n",mpi_id,all_its);
			MPI_Send(&all_its, 1, MPI_INT, 0, MSG_TAG, MPI_COMM_WORLD);
			MPI_Send(&task_its_index[0], all_its, MPI_INT, 0, MSG_TAG, MPI_COMM_WORLD);
			MPI_Send(&col_buffer[0], all_its*3, MPI_FLOAT, 0, MSG_TAG, MPI_COMM_WORLD);
//			printf("%d final send over  \n",mpi_id);
			break;
		}		
		loop++;
	}
	all_ed = rpcc();		
	time_statistics[0] = read_time; //read time	
	time_statistics[1] = athread_time;   //athread time 
	printf("read time%ld  process %ld slave time %ld, slave wait %ld\n",
			read_time/unit,  athread_time/unit, (all_ed-all_st-wait_time)/unit ,wait_time/unit);
	time_statistics[2] = wait_time;    //wait time	
	data_statistics[0] = all_photon; //photon
	data_statistics[1] = all_its; //its	
	MPI_Send(&time_statistics[0], 3, MPI_LONG, 0, MSG_TAG, MPI_COMM_WORLD);
	MPI_Send(&data_statistics[0], 2, MPI_LONG, 0, MSG_TAG, MPI_COMM_WORLD);

	free(its_buffer);
	free(task_its_index);
	free(col_list_task);
	free(cut_aabb);
	free(id_list);
}




void master_overlap(char *path, int mul){
	MPI_Status status;  
	MPI_Datatype mpi_photon;
	build_struct_photon(&mpi_photon);
	st=rpcc();
	unsigned long wait_st,wait_ed,wait_total = 0,data_dis = 0;
	int i, j, k;
	char photon_path[100];
	// data 
	read_its_data(path);
	printf("master test its num %d\n",its_num);
	float *col_buffer = (float*)malloc(its_num*3*sizeof(float));
	col_list = (float*)malloc(its_num*3*sizeof(float));	
	int *task_its_index_buffer = (int*)malloc(its_num*sizeof(int));
	printf("its num %d,",its_num);
	strcpy(photon_path,path);
	strcat(photon_path,"photonmap");  
	FILE *photon_input = fopen(photon_path, "rb");			
	fread(&m_scale, sizeof(float), 1, photon_input);
	fread(&size, sizeof(unsigned long long), 1, photon_input);
	fread(&depth, sizeof(unsigned long long), 1, photon_input);
	float aabb_box[6] = {0};
	fread(aabb_box, sizeof(float), 6, photon_input);
	fclose(photon_input);	
	printf("read over m_scale %.10f\n",m_scale);
	slave_num = numprocs-1;
	int dis_num = slave_num * mul;
	double two = 2;
	int dis_num_root = log(dis_num) / log(two);
	float child_list[6*512] = {0};
	int list_num = 0;
	printf("before bound box\n");
	get_bound_box_i(its_list, 0, its_num, aabb_box);
	printf("after bound box\n");	
	sub_bound_box(aabb_box, &aabb_box[3], child_list, &list_num, dis_num_root);

	float task_box[12]= {0};// its box + tree box
	int send_num = 0, buffer_st = 0,task_num0 = 0, subtree_num = 0;
	int max_subtree = 0,max_its = 0;
	int *power = (int*)malloc(sizeof(int) * 4) ; 
	memset(power,0,sizeof(int)*3) ; 
	ed=rpcc();
	unsigned long pre_process = (ed-st)/unit; pre_process /= 100;
	printf("pre_process %ld\n",pre_process);
	st=rpcc();
	
	for(i = 0; i < dis_num; i++){
		for (k = 0; k < 6; k++){
			task_box[k] = child_list[i * 6 + k];
		}
		float v_radius[3] = {radius,radius,radius};	
/*		if(i==3||i==11||i==7||i==13||i==6){
			v_radius[0] = 0.04179;
			v_radius[1] = 0.04179;
			v_radius[2] = 0.04179;
		}*/
		
		sub(task_box, v_radius, &task_box[6]);
		add(&task_box[3], v_radius, &task_box[9]);
		wait_st = rpcc();
		MPI_Recv(power,3,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status) ;	
//		printf("rece %d\n",power[1]);
		MPI_Send(&task_box, 12, MPI_FLOAT, power[1], MSG_TAG, MPI_COMM_WORLD);
		
		wait_ed = rpcc();
		wait_total += (wait_ed - wait_st);		
		if(power[2] == 1){
			MPI_Recv(&task_num0,1,MPI_INT,power[1],MSG_TAG,MPI_COMM_WORLD,&status) ;
			MPI_Recv(&subtree_num,1,MPI_INT,power[1],MSG_TAG,MPI_COMM_WORLD,&status) ;
			
			max_subtree = max_subtree>subtree_num?max_subtree:subtree_num;
			max_its = max_its>task_num0?max_its:task_num0;
			if(task_num0!=0){
				MPI_Recv(&task_its_index_buffer[0], task_num0, MPI_INT, power[1], MSG_TAG, MPI_COMM_WORLD,&status);
				MPI_Recv(&col_buffer[0], task_num0*3, MPI_FLOAT, power[1], MSG_TAG, MPI_COMM_WORLD,&status);
				for(j = 0;j<task_num0;j++){
					int new_col_index = task_its_index_buffer[j] * 3;
	//				printf("%d ",task_its_index_buffer[j]);
					col_list[new_col_index] = col_buffer[j*3];
					col_list[new_col_index + 1] = col_buffer[j*3+1];
					col_list[new_col_index + 2] = col_buffer[j*3+2];
				}
			}
			
		}
//		printf("rece2 %d\n",power[1]);
	}
	ed=rpcc();
	unsigned long dis = (ed-st-wait_total)/unit; 
	float task_end[12];
	task_end[0] = -1; task_end[3] = -1;  //end message
	int terminate = 0;
	st=rpcc();
	unsigned long slave_read_time = 0, slave_athread_time = 0,slave_wait_time = 0;
	printf("all send over\n");
	printf("All wait1 %ld\n  ", wait_total/unit);
	while(terminate < slave_num){
		wait_st = rpcc();
		MPI_Recv(power,3,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status) ;
		wait_ed = rpcc();
		wait_total += (wait_ed - wait_st);
		MPI_Send(&task_end, 12, MPI_FLOAT, power[1], MSG_TAG, MPI_COMM_WORLD);
		MPI_Recv(&task_num0,1,MPI_INT,power[1],MSG_TAG,MPI_COMM_WORLD,&status) ;
		MPI_Recv(&subtree_num,1,MPI_INT,power[1],MSG_TAG,MPI_COMM_WORLD,&status) ;
		max_subtree = max_subtree>subtree_num?max_subtree:subtree_num;
		max_its = max_its>task_num0?max_its:task_num0;
//		printf("its %d sutree %d\n",task_num0,subtree_num);	
		if(task_num0!=0){
			MPI_Recv(&task_its_index_buffer[0], task_num0, MPI_INT, power[1], MSG_TAG, MPI_COMM_WORLD,&status);
			MPI_Recv(&col_buffer[0], task_num0*3, MPI_FLOAT, power[1], MSG_TAG, MPI_COMM_WORLD,&status);
		}
//		printf("rece end0 %d\n",power[1]);		
		for(j = 0;j<task_num0;j++){
			int new_col_index = task_its_index_buffer[j] * 3;
//			printf("%d %d \n",j, task_its_index_buffer[j]);
			col_list[new_col_index] = col_buffer[j*3];
			col_list[new_col_index + 1] = col_buffer[j*3+1];
			col_list[new_col_index + 2] = col_buffer[j*3+2];
		}	
//		printf("rece end1 %d\n",power[1]);				
		MPI_Recv(&time_statistics[0], 3, MPI_LONG, power[1], MSG_TAG, MPI_COMM_WORLD,&status);
		slave_read_time += time_statistics[0]/100; 
		slave_athread_time += time_statistics[1]/100;
		slave_wait_time += time_statistics[2]/100;			
		terminate++;
//		printf("rece end2 %d\n",power[1]);
		
		
	}
	ed=rpcc();
	dis += (ed-st-wait_total)/unit; 
	dis /= 100;
	printf("max its %d max sub tree %d\n", max_its, max_subtree);
	printf("distribute %ld\n",dis);
	printf("All wait2 %ld\n  ", wait_total/unit);
	

	char reslut_path[100], str[25];
	strcpy(reslut_path, path);
	strcat(reslut_path,"timestatic");	
	FILE *resoutput = fopen(reslut_path,"a");	
	fprintf(resoutput, "max its %d max sub tree %d distribute %ld  wait %ld \n",max_its, max_subtree,dis,wait_total/unit);
	
	
	
	printf("slave read %ld athread %ld wait %ld \n",slave_read_time,slave_athread_time,slave_wait_time);
}


void slave_overlap(char *path){
	MPI_Status status;  
	MPI_Datatype mpi_photon;
	build_struct_photon(&mpi_photon);
	unsigned long read_time = 0, athread_time = 0,wait_time = 0, read_st, read_ed, ath_st, ath_ed;
	int i, j, k, slave_photon_size;
	int *power ; 
	power = (int*)malloc(sizeof(int)*3) ;  
	power[0] = 0 ;  power[1] = mpi_id ; 
	id_list = (int*)malloc(1*sizeof(int));
	createPrecompTables();	
	//read all its
	read_its_data(path);
	task_its_list=(struct Intersection*)malloc(its_num*sizeof(struct Intersection));
	task_its_index = (int*)malloc(its_num*sizeof(int));	
	col_list_task=(float*)malloc(its_num*4*sizeof(float));	
	
	///read photon option data 	
	st = rpcc();
	char photon_path[100];
	strcpy(photon_path,path);
	strcat(photon_path,"photonmap");  
	FILE *photon_input = fopen(photon_path, "rb");
	fread(&m_scale, sizeof(float), 1, photon_input);
	fread(&size, sizeof(unsigned long long), 1, photon_input);
	fread(&depth, sizeof(unsigned long long), 1, photon_input);
	int max_subtree = 200000000;   //5.21G
	sub_tree=(struct Photon*)malloc(max_subtree * sizeof(struct Photon));
	float aabb_box[6] = {0};
	fread(aabb_box, sizeof(float), 6, photon_input);
	
	float task_box[12] = {0};  //its box  photon box
	MPI_Send(power, 3, MPI_INT, 0, MSG_TAG,MPI_COMM_WORLD);
	MPI_Recv(&task_box, 12, MPI_FLOAT, 0, MSG_TAG, MPI_COMM_WORLD, &status);
/*	printf("mpi %d \n",mpi_id);
	for(i = 0;i<12;i++){
		printf("%f ",task_box[i]);
	}
	printf("\n");  */
	int seek_step = 15;
	for(;;){
		//read photon 
		float pos_temp[3];
		int task_its_num = 0, subtree_num = 0;
		st = rpcc(); 
		for (i = 0; i < size; i++){ 		
			fread(&sub_tree[subtree_num].pos[0], sizeof(float), 3, photon_input);
			if(in_sub_tree(&task_box[6], sub_tree[subtree_num].pos)==1){
				//read
				fread(&sub_tree[subtree_num].right, sizeof(unsigned int), 1, photon_input);
				fread(&sub_tree[subtree_num].data[0], sizeof(unsigned char), 8, photon_input);
				fread(&sub_tree[subtree_num].dep, sizeof(unsigned short), 1, photon_input);
				fread(&sub_tree[subtree_num].flag, sizeof(char), 1, photon_input);
				if(subtree_num>=max_subtree-1){
					printf("sub tree too big \n");
					break;
				}
				subtree_num++;
				
			}else{
				fseek(photon_input, seek_step, SEEK_CUR);
			}
		}
		for (i = 0; i < its_num; i++){ 		
			if (in_sub_tree(task_box, its_list[i].p)==1){
				task_its_list[task_its_num] = its_list[i];
				task_its_index[task_its_num] = i;
				task_its_num++;
			}
		}
		construct_photon_tree_xx(1, sub_tree, 0, subtree_num - 1, &task_box[9], &task_box[6]);
//		printf("subtree_num %d %d task_its_num %d \n",subtree_num, power[1],task_its_num);
		task_num[0] = task_its_num;
		task_num[1] = subtree_num;
		task_num[2] = 1;
		id_list[0] = 0;
		int type = 1;
		cut_aabb = task_box;	
		ed = rpcc(); 
		read_time += ed-st;
		st = rpcc(); 
		athread_spawn(NNsearch,&type);
		athread_join();
		ed = rpcc(); 
		athread_time += ed-st;
		power[2]=1;
		st = rpcc(); 
		MPI_Send(power, 3, MPI_INT, 0, MSG_TAG,MPI_COMM_WORLD) ;
		MPI_Recv(&task_box, 12, MPI_FLOAT, 0, MSG_TAG, MPI_COMM_WORLD, &status);
		MPI_Send(&task_its_num, 1, MPI_INT, 0, MSG_TAG, MPI_COMM_WORLD);
		MPI_Send(&subtree_num, 1, MPI_INT, 0, MSG_TAG, MPI_COMM_WORLD);
		ed = rpcc(); wait_time += ed-st;		
		if(task_its_num!=0){	
//			printf("send col\n");
			MPI_Send(&task_its_index[0], task_num[0], MPI_INT, 0, MSG_TAG, MPI_COMM_WORLD);
			MPI_Send(&col_list_task[0], task_its_num*3, MPI_FLOAT, 0, MSG_TAG, MPI_COMM_WORLD);		
		}				
		if(task_box[0]==-1&&task_box[3] == -1){
			time_statistics[0] = read_time/unit; //read time	
			time_statistics[1] = athread_time/unit;   //athread time 
			time_statistics[2] = wait_time/unit;    //estimate time	
			MPI_Send(&time_statistics[0], 3, MPI_LONG, 0, MSG_TAG, MPI_COMM_WORLD);
			break;
		}	
		
	}
	fclose(photon_input);
	
	free(task_its_list);
	free(task_its_index);
	free(col_list_task);
}


void master_disk_Kato(char *path){
	int k,i,j;	
	char photon_path[100];
	strcpy(photon_path,path);
	strcat(photon_path,"photonmap");  
	FILE *photon_input = fopen(photon_path, "rb");		
	fread(&m_scale, sizeof(float), 1, photon_input);
	fread(&size, sizeof(unsigned long long), 1, photon_input);
	fclose(photon_input);
	
	read_its_data(path);
	float *col_buffer = (float*)malloc(its_num*4*sizeof(float));
	col_list = (float*)malloc(its_num*3*sizeof(float));	
	float *radius_list = (float*)malloc(its_num*sizeof(float));
	for(i = 0;i<its_num;i++){
		col_list[i*3] = 0;col_list[i*3+1] = 0;col_list[i*3+2] = 0;
		col_buffer[i*4] = 0;col_buffer[i*4+1] = 0;col_buffer[i*4+2] = 0;col_buffer[i*4+3] = 0;
		radius_list[i] = 0;
	}
	
	MPI_Status status ; 
	MPI_Status Istatus[10] ; 
	MPI_Datatype mpi_photon;
	build_struct_photon(&mpi_photon);
	unsigned long wait_st,wait_ed,wait_total = 0, data_dis = 0, all_st, all_ed;
	double statistics[5][256], time_static[10];	
	long long all_photon = 0;	
	
	st=rpcc();
	int *power = (int*)malloc(sizeof(int) * 3) ; 
	memset(power,0,sizeof(int)*3) ; 
	all_st = rpcc();	

	int terminate = 0;
	wait_st = rpcc();
	slave_num = numprocs-1;
	while(terminate < slave_num){		
		MPI_Recv(power, 3, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status) ;
		printf("power %d \n",power[1]);
		MPI_Recv(&col_buffer[0], its_num*4, MPI_FLOAT, power[1], MSG_TAG, MPI_COMM_WORLD,&status);
		wait_ed = rpcc();	
		wait_total += (wait_ed - wait_st);	
		for(j = 0;j<its_num;j++){
			int new_col_index = j * 3;
			col_list[j * 3] += col_buffer[j*4];
			col_list[j * 3 + 1] += col_buffer[j*4+1];
			col_list[j * 3 + 2] += col_buffer[j*4+2];
			radius_list[j] += col_buffer[j*4+3];
		}
		data_dis += (rpcc() - wait_st);	
		
		MPI_Recv(&time_statistics[0], 3, MPI_LONG, power[1], MSG_TAG, MPI_COMM_WORLD,&status);
		
		statistics[0][terminate] = time_statistics[0]/100; 
		statistics[1][terminate] = time_statistics[1]/100; 
		statistics[2][terminate] = time_statistics[2]/100;		
		statistics[3][terminate] = data_statistics[0]; 
		statistics[4][terminate] = data_statistics[1];				
		terminate++;
		wait_st = rpcc();
	}
	for(j = 0;j<its_num;j++){
		dot_vf(&col_list[j * 3], slave_num * m_scale * 3 * INV_PI / radius_list[j], &col_list[j * 3]);
	}	
	ed=rpcc();
	time_static[0] = (ed-st)/unit;time_static[0]/=100;				

	st=rpcc();
	float task_box[6]= {0};
	int task_num0 = 0, subtree_num = 0;
	ed=rpcc();
	time_static[1] = (ed-st-wait_total-data_dis)/unit; time_static[1] /= 100;	
	st=rpcc();
	

	all_ed = rpcc();
	time_static[2] = data_dis/unit; time_static[2] /= 100;
	time_static[3] = wait_total/unit; time_static[3] /= 100;
	time_static[4] = (all_ed-all_st)/unit; time_static[4] /= 100; 
	print_static(statistics, time_static,all_photon, numprocs, 1, path, 0);
	printf("All wait %ld\n  ", wait_total/unit);
	printf("All Communication %ld\n  ", data_dis/unit);
	free(power);
}

void slave_disk_Kato(char *path){   //read from disk  send to master  
	MPI_Status status;  
	MPI_Datatype mpi_photon;
	build_struct_photon(&mpi_photon);
	///read photon option data 	
	int i, j, k, slave_photon_size;
	st = rpcc();
	char photon_path[100], its_path[100];
	strcpy(photon_path,path);
	strcat(photon_path,"photonmap");  
	FILE *photon_input = fopen(photon_path, "rb");		
	
	fread(&m_scale, sizeof(float), 1, photon_input);
	fread(&size, sizeof(unsigned long long), 1, photon_input);
	fread(&depth, sizeof(unsigned long long), 1, photon_input);
	float aabb_box[6] = {0};
	fread(aabb_box, sizeof(float), 6, photon_input);
			printf("%d read ok",mpi_id);
	slave_num = numprocs-1;
	slave_photon_size = size/slave_num;	
	sub_tree=(struct Photon*)malloc(slave_photon_size*sizeof(struct Photon));
	id_list = (int*)malloc(1*sizeof(int));
	fseek(photon_input, (mpi_id-1) * (sizeof(struct Photon)-1), SEEK_CUR);
	long curpos, length;
	int step = (slave_num-1) * (sizeof(struct Photon)-1);
	for (i = 0; i < slave_photon_size; i++){ 		
		fread(&sub_tree[i].pos[0], sizeof(float), 3, photon_input);
		fread(&sub_tree[i].right, sizeof(unsigned int), 1, photon_input);
		fread(&sub_tree[i].data[0], sizeof(unsigned char), 8, photon_input);
		fread(&sub_tree[i].dep, sizeof(unsigned short), 1, photon_input);
		fread(&sub_tree[i].flag, sizeof(char), 1, photon_input);
		fseek(photon_input, step, SEEK_CUR);
		curpos=ftell(photon_input);	
	}
	fclose(photon_input);
	
	read_its_data(path);
	task_its_list = its_list;	
	col_list_task=(float*)malloc(its_num*4*sizeof(float));	
	createPrecompTables();	
	ed = rpcc();
	time_statistics[0] = (ed - st)/unit; //read disk time 
//	printf("%d read ok3",mpi_id);
	float *aabb = (float*)malloc(6*sizeof(float));
	st = rpcc();
	get_bound_box_p(sub_tree, 0, slave_photon_size, aabb);	
	construct_photon_tree_xx(1, sub_tree, 0, slave_photon_size - 1, &aabb[3], aabb);
	ed = rpcc();
	time_statistics[1] = (ed - st)/unit;   //tree time 
	
	st = rpcc();
	task_num[0] = its_num;
	task_num[2] = 1;
	id_list[0] = 0;
	int type = 2;
	cut_aabb = aabb;	
	athread_spawn(NNsearch,&type);
	athread_join();
	ed = rpcc();	
	time_statistics[2] =(ed-st)/unit;
	////send data
	int *power ; 
	power = (int*)malloc(sizeof(int)*3) ;  
	power[0] = 0 ;  power[1] = mpi_id ;  // [0] end  [1] id  [2] if have work 
	MPI_Send(power, 3, MPI_INT, 0, MSG_TAG,MPI_COMM_WORLD) ;
	MPI_Send(&col_list_task[0], its_num*4, MPI_FLOAT, 0, MSG_TAG, MPI_COMM_WORLD);
	//end 	
	MPI_Send(&time_statistics[0], 3, MPI_LONG, 0, MSG_TAG, MPI_COMM_WORLD);
	
	free(task_its_list);
	free(sub_tree);
	free(col_list_task);
	return;
}
void only_one(char *path, float *aabb_box){
	MPI_Status status;
	MPI_Datatype mpi_photon;
	build_struct_photon(&mpi_photon);
	int i = 0,all_index = 0, list_num = 0;
	int all_photon = 0;
	
	unsigned long st1, ed1, distribute_time = 0, wait_time = 0 , true_distribute = 0;	
	printf("start construct\n");
	st=rpcc(); 
	construct_photon_tree_xx(1, x_tree, 0, size, &aabb_box[3], aabb_box);
	ed=rpcc();
	printf("construct photon tree %ld\n", (ed-st)/unit);
	st=rpcc(); 
	int type = 1;
	task_num[0] = its_num;

	free(sub_tree);	
	free(task_its_list);
	free(task_its_index);
	free(col_list_task);
		
	sub_tree = x_tree;
	task_its_list = its_list;
	col_list_task = col_list;

	athread_spawn(NNsearch,&type);
	athread_join();
	printf("%d manycore counter=%ld subtree size %d task num %d\n", mpi_id, ed-st, task_num[1], task_num[0]);
	ed=rpcc();
	printf("estimateIrradiance time %ld\n", (ed-st)/unit);
	free(x_tree);
	free(col_list);	
	free(its_list);
}
int CharToInt(char *a){
	if(a[0]=='1'){
		if(a[1]=='6')
			return 16;
		else if(a[1]=='2')
			return 128;
		else
			return 1;
	}
	if(a[0]=='2'){
		return 2;
	}
	if(a[0]=='3'){
		if(a[1]=='2')
			return 32;
		return 3;
	}
	if(a[0]=='4')
		return 4;
	if(a[0]=='6')
		return 64;
	if(a[0]=='8'){
		return 8;
	}
	printf("char to int error\n");
	return 1;
}
struct Morton_node * gen_morton_code(char *path){
	char photon_path[100];
	strcpy(photon_path,path);
	strcat(photon_path,"photonmap");  
	FILE *photon_input = fopen(photon_path, "rb");		
	float aabb[6] = {0};
	
	fread(&m_scale, sizeof(float), 1, photon_input);
	fread(&size, sizeof(unsigned long long), 1, photon_input);
	fread(&depth, sizeof(unsigned long long), 1, photon_input);	
	fread(aabb, sizeof(float), 6, photon_input);
	
	int read_length = size/numprocs;
	fseek(photon_input, mpi_id * read_length * (sizeof(struct Photon)-1), SEEK_CUR);
	sub_tree=(struct Photon*)malloc((read_length+10) *sizeof(struct Photon));
	struct Morton_node *morton_list;
//	if(mpi_id==0)
//		morton_list = (struct Morton_node*)malloc((size+10) *sizeof(struct Morton_node));	
//	else
		morton_list = (struct Morton_node*)malloc((read_length+10) *sizeof(struct Morton_node));
	long curpos;	
	float aabb_range[3] = { aabb[3] - aabb[0], aabb[4] - aabb[1], aabb[5] - aabb[2]};
	int i, base_index = mpi_id * read_length;
  	for (i = 0; i < read_length; i++){ 	

		// 这里可以只读pos就行不用全部读取
	
		fread(&sub_tree[i].pos[0], sizeof(float), 3, photon_input);		
		fread(&sub_tree[i].right, sizeof(unsigned int), 1, photon_input);
		fread(&sub_tree[i].data[0], sizeof(unsigned char), 8, photon_input);
		fread(&sub_tree[i].dep, sizeof(unsigned short), 1, photon_input);
		fread(&sub_tree[i].flag, sizeof(char), 1, photon_input);
		/////这里的i 是相对的i 
		morton_list[i].its_id = base_index + i;
//		morton_list[i].code = mortonid(sub_tree[i], aabb, aabb_range);
		morton_list[i].code = mortonid_ph(sub_tree[i], aabb, aabb_range);
	}
	fclose(photon_input);
	return morton_list;
}



void mpi_psrs_phase1(struct Morton_node *morton_list, int sort_task_num, struct Morton_node *main_e){
	struct Morton_node *main_e_rec = (struct Morton_node*)malloc(numprocs*numprocs*sizeof(struct Morton_node));
	int i;
	quickSort(morton_list, 0, sort_task_num - 1);

	//subtree也重新排序
	struct Photon *new_subtree = (struct Photon*)malloc((sort_task_num+10) *sizeof(struct Photon));
	int base_id = mpi_id * sort_task_num;
	for (i = 0; i < sort_task_num; i++){
		new_subtree[i] = sub_tree[morton_list[i].its_id - base_id];
		morton_list[i].its_id = base_id + i;
	}
	free(sub_tree);
	sub_tree = new_subtree;
		
	//select main element	
	int sample_num = sort_task_num/numprocs;
	for (i = 0; i < numprocs; i++) {
		main_e[i].code = morton_list[i * sample_num].code;    
	}
	//gather main element and sort 
	MPI_Gather(main_e,numprocs*2,MPI_INT,main_e_rec,numprocs*2,MPI_INT,0,comm_psrs); 
	if(mpi_id==0){
		quickSort(main_e_rec, 0, numprocs*numprocs - 1);
		for(i = 0;i<numprocs-1;i++){
			main_e[i] = main_e_rec[numprocs*(i+1)];
//			printf("%d  ",main_e[i].code);
		}
//		printf("\n");
		main_e[numprocs-1].code = INT_MAX;
	}
	MPI_Bcast(main_e, 2*numprocs, MPI_INT, 0, comm_psrs);
	free(main_e_rec);
}
struct BVH_tree * mpi_psrs_phase2(struct Morton_node *morton_list, int sort_task_num, 
					struct Morton_node *main_e, char *path){
	MPI_Datatype mpi_photon,mpi_bvh_tree;
	build_struct_photon(&mpi_photon);
	
	
	int *partitionSizes = (int*)malloc(numprocs*sizeof(int));
	int i,j;
	for ( i = 0; i < numprocs; i++) {
		partitionSizes[i] = 0;
    }
	int index = 0;
	//根据主元计算划分
	for ( i = 0; i < sort_task_num; i++) {
		while (morton_list[i].code > main_e[index].code) {
			index += 1;
		}
		if (index == numprocs-1) {
		  partitionSizes[numprocs - 1] = sort_task_num - i;
		  break;
		}
		partitionSizes[index]++ ;   //划分大小自增
	}
	int *newPartitionSizes = (int*)malloc(numprocs*sizeof(int));	
	int totalSize = 0;
	int *sendDisp = (int *) malloc(numprocs * sizeof(int));
	int *recvDisp = (int *) malloc(numprocs * sizeof(int));	
	MPI_Alltoall(partitionSizes, 1, MPI_INT, newPartitionSizes, 1, MPI_INT, comm_psrs);		
	// 计算划分的总大小，并给新划分分配空间
	for ( i = 0; i < numprocs; i++) {
		totalSize += newPartitionSizes[i];
	}
	sendDisp[0] = 0;
	recvDisp[0] = 0;
	for ( i = 1; i < numprocs; i++) {
		sendDisp[i] = partitionSizes[i - 1] + sendDisp[i - 1];
		recvDisp[i] = newPartitionSizes[i - 1] + recvDisp[i - 1];
	}
	////sub tree alltoall 传输
//	printf("%d  before malloc \n",mpi_id);
	struct Photon *new_subtree = (struct Photon*)malloc((totalSize+10) *sizeof(struct Photon));
//	printf("%d before allotall \n",mpi_id);
	MPI_Alltoallv(&sub_tree[0], partitionSizes, sendDisp, mpi_photon, new_subtree, newPartitionSizes, recvDisp, mpi_photon, comm_psrs);
	free(sub_tree);
	sub_tree = new_subtree;
//	printf("%d after allotall \n",mpi_id);
	for(i = 0;i<numprocs;i++){
		sendDisp[i]*=2;
		recvDisp[i]*=2;
		partitionSizes[i]*=2;
		newPartitionSizes[i]*=2;
	}
	
	struct Morton_node *newPartitions = (struct Morton_node *) malloc(totalSize * sizeof(struct Morton_node));
	MPI_Alltoallv(&morton_list[0], partitionSizes, sendDisp, MPI_INT, newPartitions, newPartitionSizes, recvDisp, MPI_INT, comm_psrs);
//	printf("%d after Morton_node allotall \n",mpi_id);
	for(i = 0;i<numprocs;i++){
		recvDisp[i]/=2;
		newPartitionSizes[i]/=2;
	}	
	int *indexes, *partitionEnds, *subListSizes, totalListSize;
	indexes = (int *) malloc(numprocs * sizeof(int));
	partitionEnds = (int *) malloc(numprocs * sizeof(int));
	indexes[0] = 0;
	totalListSize = newPartitionSizes[0];
	for ( i = 1; i < numprocs; i++) {
		totalListSize += newPartitionSizes[i];
		indexes[i] = indexes[i-1] + newPartitionSizes[i-1];
		partitionEnds[i-1] = indexes[i];
	}
	partitionEnds[numprocs - 1] = totalListSize;
	struct Morton_node *sortedSubList = (struct Morton_node *) malloc(totalListSize * sizeof(struct Morton_node));
	subListSizes = (int *) malloc(numprocs * sizeof(int));
	
	//生成新的文件 命名为 photon_1、photon_2
	// slave直接建树  kd树 计算时可以直接拿来用
	float *aabb = (float *) malloc(6 * sizeof(float));
//	float *aabb_list = (float *) malloc(6 * numprocs * sizeof(float));
	get_bound_box_p(sub_tree, 0, totalListSize, aabb);
	construct_photon_tree_xx(1, sub_tree, 0, totalListSize - 1, &aabb[3], aabb);	
//	printf("%d %f %f %f %f %f %f\n", mpi_id, aabb[0],aabb[1],aabb[2],aabb[3],aabb[4],aabb[5]);
	char new_path[100], photon_id[10];
	strcpy(new_path,path); strcat(new_path,"new_photon_"); 
	sprintf(photon_id, "%d", mpi_id);
	strcat(new_path,photon_id); 	
	FILE *photon_new = fopen(new_path, "wb");	
//	fwrite(&sub_tree, sizeof(struct Photon), totalListSize, photon_new);
	for ( i = 0; i < totalListSize; i++) {
		if(fwrite(&sub_tree[i], sizeof(struct Photon), 1, photon_new)!=1)
             printf("file write error\n");
	}	
//	printf("totalListSize  %d \n",totalListSize);
	fclose(photon_new);	
	
	// slave把 自己的bbox以及 部分上层节点发送给master 
//	for(i = 0;i<6;i++){
//		printf("%f ",aabb[i]);
//	}
//	printf("\n");
//	MPI_Gather(aabb, 6, MPI_FLOAT, aabb_list, 6, MPI_FLOAT, 0, comm_psrs);
	
	double two = 2;
	int numprocs_root = log(numprocs) / log(two);
	int photon_tree_lv = BVH_TREE_LV;      //top tree 总层数
	int slave_top_tree_lv = photon_tree_lv - numprocs_root;   //    每个slave的toptree的层数
	int slave_top_tree_size = pow(2.0, slave_top_tree_lv) - 1;      // slave top tree size 
	
	
	struct BVH_tree *top_tree = (struct BVH_tree*)malloc(slave_top_tree_size *sizeof(struct BVH_tree));	
	struct BVH_tree *top_tree_list = (struct BVH_tree*)malloc(numprocs * (slave_top_tree_size+3) *sizeof(struct BVH_tree));
	
	get_top_tree(sub_tree, top_tree, slave_top_tree_lv, aabb, mpi_id, slave_top_tree_size, totalListSize);
	build_struct_bvh(&mpi_bvh_tree, top_tree);
	
		
	MPI_Gather(top_tree, slave_top_tree_size, mpi_bvh_tree, top_tree_list, slave_top_tree_size, mpi_bvh_tree, 0, comm_psrs);
	MPI_Gather(&totalListSize, 1, MPI_INT, subListSizes, 1, MPI_INT, 0, comm_psrs);	

	if (mpi_id == 0) {
		construct_photon_bvh(subListSizes, top_tree_list, numprocs, slave_top_tree_size);			
	}
	return top_tree_list;
}
struct BVH_tree * mpi_psrs(char *path){
	MPI_Status status ; 
	struct Morton_node *morton_list = gen_morton_code(path);
//	printf("mpi %d gen code ok\n",mpi_id);
	int sort_task_num = size / numprocs;
	int sample_num = sort_task_num / numprocs;
	struct Morton_node *main_e = (struct Morton_node*)malloc(numprocs*sizeof(struct Morton_node)); 
	mpi_psrs_phase1(morton_list, sort_task_num, main_e);	
//	printf("mpi %d  psrs phase1 ok\n",mpi_id);
	struct BVH_tree *photon_bvh_tree = mpi_psrs_phase2(morton_list, sort_task_num, main_e, path);
//	printf("mpi %d  psrs phase2 ok\n",mpi_id);
	return photon_bvh_tree;
}


int mpi_psrs_its(char *path, int mul, struct Morton_node *morton_list,
				int *its_queue_range,float *its_queue_aabb ){
					
	MPI_Status status ; 
	read_its_data(path);
	float aabb[6] = {0};
	get_bound_box_i(its_list, 0, its_num, aabb);
	float aabb_range[3] = { aabb[3] - aabb[0], aabb[4] - aabb[1], aabb[5] - aabb[2]};
	int i,j,read_num = 0;
	
/*	for(i = 0;i<6;i++) {
		printf("%f ", aabb[i]);
	}
	printf("\n");
	printf("%f %f %f\n",aabb_range[0],aabb_range[1],aabb_range[2]);*/
	
	int read_length = its_num/numprocs;
	int st = read_length*mpi_id;
//	printf("%d %d\n", read_length, st);
	
	for (i = 0; i < read_length; i++){ 	
		int temp_code = mortonid_its(its_list[st+i], aabb, aabb_range);
		if(temp_code!=0){
			morton_list[read_num].its_id = i+st;
			morton_list[read_num].code = temp_code;
			read_num++;
		}
	}
	read_length = read_num;
//	printf("read_length %d\n",read_length);
//	printf("code ok\n");
	quickSort(morton_list, 0, read_length - 1);
//	printf("sort ok\n");
	struct Morton_node *main_e = (struct Morton_node*)malloc(numprocs*sizeof(struct Morton_node)); 	
	struct Morton_node *main_e_rec = (struct Morton_node*)malloc(numprocs*numprocs*sizeof(struct Morton_node));
	
	int sample_num = read_length / numprocs;
	for (i = 0; i < numprocs; i++) {
		main_e[i].code = morton_list[i * sample_num].code; 
//			printf("%d ",main_e[i].code);
	}
//	printf("\n");
//	printf("slave main e ok\n");
	//gather main element and sort
	int sdrv_main = numprocs*2;
	MPI_Gather(main_e, sdrv_main ,MPI_INT,main_e_rec, sdrv_main ,MPI_INT,0,MPI_COMM_WORLD); 
	if(mpi_id==0){
		quickSort(main_e_rec, 0, numprocs * numprocs - 1);
		for(i = 0;i<numprocs-1;i++){
			main_e[i] = main_e_rec[numprocs*(i+1)];
//			printf("%d ",main_e[i].code);
		}
//		printf("\n\n\n");
		main_e[numprocs-1].code = INT_MAX;
	}
	MPI_Bcast(main_e, 2*numprocs, MPI_INT, 0, MPI_COMM_WORLD);
	
	int *partitionSizes = (int*)malloc(numprocs*sizeof(int));
	for ( i = 0; i < numprocs; i++) {
		partitionSizes[i] = 0;
    }
	int index = 0;
	//根据主元计算划分
	for ( i = 0; i < read_length; i++) {
		while (morton_list[i].code > main_e[index].code) {
			index += 1;
		}
		if (index == numprocs-1) {
		  partitionSizes[numprocs - 1] = read_length - i;
		  break;
		}
		partitionSizes[index]++ ;   //划分大小自增
	}
	
	int *newPartitionSizes = (int*)malloc(numprocs*sizeof(int));	
	int totalSize = 0;
	int *sendDisp = (int *) malloc(numprocs * sizeof(int));
	int *recvDisp = (int *) malloc(numprocs * sizeof(int));	
	MPI_Alltoall(partitionSizes, 1, MPI_INT, newPartitionSizes, 1, MPI_INT, MPI_COMM_WORLD);		
	// 计算划分的总大小，并给新划分分配空间
	for ( i = 0; i < numprocs; i++) {
		totalSize += newPartitionSizes[i];
	}
//	printf("%d %d \n", mpi_id,totalSize);
	sendDisp[0] = 0;
	recvDisp[0] = 0;
	for ( i = 1; i < numprocs; i++) {
		sendDisp[i] = partitionSizes[i - 1] + sendDisp[i - 1];
		recvDisp[i] = newPartitionSizes[i - 1] + recvDisp[i - 1];
	}
	for ( i = 0; i < numprocs; i++) {
		sendDisp[i] *= 2;
		recvDisp[i] *= 2;
		partitionSizes[i] *= 2;
		newPartitionSizes[i] *= 2;
	}
	// morton code all to all
	struct Morton_node *new_morton_list = (struct Morton_node*)malloc((totalSize+10) * sizeof(struct Morton_node));
	MPI_Alltoallv(&morton_list[0], partitionSizes, sendDisp, MPI_INT, new_morton_list, 
							newPartitionSizes, recvDisp, MPI_INT, MPI_COMM_WORLD);
	int *SlaveSizesList = (int *) malloc(numprocs * sizeof(int));
//	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allgather(&totalSize, 1, MPI_INT, SlaveSizesList, 1, MPI_INT, MPI_COMM_WORLD);
//	printf("MPI_Allgather \n");
	// slave sort  可能需要改成多路归并
	quickSort(new_morton_list, 0, totalSize - 1);
	
/*	if(mpi_id==1){
		for(i = 0;i<totalSize;i++)
			printf("%d ",new_morton_list[i].code);
		printf("\n");
	}*/
//	printf("sort over \n");
	// master gather all morton code
	// 计算各个进程上的相对于recvbuf的偏移量	
	recvDisp[0] = 0;
	int valid_its_num = 0;
	for ( i = 0; i < numprocs; i++) {
		valid_its_num+=SlaveSizesList[i];
		SlaveSizesList[i] *= 2;
//		printf("SlaveSizesList[i] %d\n",SlaveSizesList[i]);
		if(i!=0)
			recvDisp[i] = SlaveSizesList[i - 1] + recvDisp[i - 1];
//		printf("SrecvDisp[i]] %d\n",recvDisp[i]);

	}
//	MPI_Barrier(MPI_COMM_WORLD);	
	MPI_Allgatherv(&new_morton_list[0], totalSize*2, MPI_INT, morton_list, SlaveSizesList, recvDisp, MPI_INT, MPI_COMM_WORLD);
	int its_queue_num = 0;
//	printf("MPI_Allgatherv\n");
	if(mpi_id==0){
		int max = 1024*1024*1024; //new_morton_list[totalSize-1].code
		morton_list[valid_its_num].code = max;
		int min = morton_list[0].code;
	
		
		int change_point = (numprocs-1)*mul;// /2;
		int head = 0, rear = valid_its_num/(numprocs-1)/mul; //32769;  //32768;//
		printf("max its_num %d\n",rear);
		construct_morton_queue(its_list, morton_list, 
						min, max, its_queue_aabb, its_queue_range,
						valid_its_num, &its_queue_num, &head, &rear, &change_point);
	
/*		for(i = 0;i<its_queue_num;i++){
					printf("%d  ", its_queue_range[i]);
					int k;
					for (k = 0;k<6;k++)
						printf("  %f  ",its_queue_aabb[i*6 + k]);
					printf("\n");
		}*/
	}
	
	free(new_morton_list);

//	printf("morton list %d  %d\n",morton_list[0].code,morton_list[0].its_id);
	return its_queue_num;
}

int main(int argc, char *argv[])
{
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&mpi_id);
	
	slave_num = numprocs-1;
#ifdef USE_SLAVE
	slave_num = USE_SLAVE;
#endif
	
	MPI_Group world_group;
	MPI_Comm_group(MPI_COMM_WORLD, &world_group);
	
	
	int ranks[1] = {numprocs-1};	
	MPI_Group_excl(world_group, 1, ranks, &G_psrs);
	MPI_Comm_create(MPI_COMM_WORLD, G_psrs, &comm_psrs);
	
	
	athread_init();
	unit = 15000000;//0.01s
	task_num[2] = 1;
	unsigned long tree_st,tree_ed;
	
	
	int ph_bvh_size = pow(2.0, BVH_TREE_LV) - 1; 
	struct BVH_tree *photon_bvh_tree = NULL;
	struct Morton_node *morton_list = (struct Morton_node*)malloc(1024*1024*32 * sizeof(struct Morton_node));
	int *its_queue_range = (int*)malloc(6000*2*sizeof(int));
	float *its_queue_aabb = (float*)malloc(6000*6*sizeof(float));

	
	if(numprocs==1){  //only one node master
		printf("start\n");
		float aabb_box[6] = {0};
		read_photon_data(argv[1], aabb_box);
		read_its_data(argv[1]);
		printf("read over\n");
		all_st=rpcc();		
		only_one(argv[1], aabb_box);				
		all_ed = rpcc();
		printf("only one total time is %ld\n", (all_ed - all_st)/unit);
		write_color(argv[1]);
		free_all();
	}
	else{
		int mul = CharToInt(argv[3]);  //task mutiple
		if(mpi_id==0){
			printf("start\n");				
			all_st=rpcc();		
			if(!strcmp(argv[2], "C")){  	//kato 2001					
				master_disk_Kato(argv[1]);	
			}
			else if(!strcmp(argv[2], "O")){  //overlap			
				master_overlap(argv[1], mul);
			}
			else { 	
				printf("mul %d \n",mul);			
				tree_st=rpcc();
				numprocs -=1;
				photon_bvh_tree = mpi_psrs(argv[1]);	
				MPI_Datatype mpi_bvh_tree;
				build_struct_bvh(&mpi_bvh_tree, photon_bvh_tree);
				MPI_Bcast(photon_bvh_tree, ph_bvh_size, mpi_bvh_tree, 0, MPI_COMM_WORLD);
//				printf("test bcast  %f   %d %d",photon_bvh_tree[2].box[1],photon_bvh_tree[2].children[1],photon_bvh_tree[2].range[1]);

				tree_ed=rpcc();
				printf("tree time is %ld\n", (tree_ed - tree_st)/unit);
				numprocs +=1;
							
				tree_st=rpcc();
				int its_task_num = mpi_psrs_its(argv[1], mul, morton_list, its_queue_range, its_queue_aabb);
				printf("its_task_num %d \n",its_task_num);
				int i;
/*				for(i = 0;i<its_task_num;i++){
					printf("%d  ", its_queue_range[i]);
					int k;
					for (k = 0;k<6;k++)
						printf("  %f  ",its_queue_aabb[i*6 + k]);
					printf("\n");
				}*/
				
//				int task_queue[520],rec_queue[520], displs[520],range_queue[520];
				float *queue_aabb = (float*)malloc(1024*6*sizeof(float));
		
				tree_ed=rpcc();
				printf("its tree time is %ld\n", (tree_ed - tree_st)/unit);				
				tree_st=rpcc();
				master_pre_estimation(morton_list, its_queue_aabb, its_queue_range, its_task_num, photon_bvh_tree, argv[1]);
//				master_pre_estimation(morton_list, queue_aabb,range_queue, task_num, photon_bvh_tree,argv[1]);

			}
			tree_st=rpcc();
//			normalize_color();
			write_color(argv[1]);
			tree_ed=rpcc();
			printf("wirte time is %ld\n", (tree_ed - tree_st)/unit);
//			free_all();
			all_ed = rpcc();
			printf("total time is %ld\n", (all_ed - all_st)/unit);
			char reslut_path[100], str[25];
			strcpy(reslut_path, argv[1]);
			strcat(reslut_path,"timestatic");	
			FILE *resoutput = fopen(reslut_path,"a");	
			fprintf(resoutput, "total  %ld\n",(all_ed - all_st)/unit);
		}
		else{						
			if(!strcmp(argv[2], "C")){
				//slave(2);				
				slave_disk_Kato(argv[1]);
			}
			else if(!strcmp(argv[2], "O")){
				slave_overlap(argv[1]);
			}
			else{
				
				numprocs -=1;
				if(mpi_id!= numprocs){
					mpi_psrs(argv[1]);  //			
				}					
				else{
					char photon_path[100];
					strcpy(photon_path,argv[1]);
					strcat(photon_path,"photonmap");  
					FILE *photon_input = fopen(photon_path, "rb");						
					fread(&m_scale, sizeof(float), 1, photon_input);
					fread(&size, sizeof(unsigned long long), 1, photon_input);
					fread(&depth, sizeof(unsigned long long), 1, photon_input);	
					fclose(photon_input);
				}
				numprocs +=1;
				MPI_Datatype mpi_bvh_tree;
				build_struct_bvh(&mpi_bvh_tree, photon_bvh_tree);
				photon_bvh_tree = (struct BVH_tree *)malloc(ph_bvh_size*sizeof(struct BVH_tree)) ;
				MPI_Bcast(photon_bvh_tree, ph_bvh_size, mpi_bvh_tree, 0, MPI_COMM_WORLD);
//				printf("test bcast  %f   %d %d",photon_bvh_tree[2].box[1],photon_bvh_tree[2].children[1],photon_bvh_tree[2].range[1]);
				
				int its_task_num = mpi_psrs_its(argv[1], mul, morton_list, its_queue_range, its_queue_aabb);
//				printf("morton list %d  %d\n",morton_list[0].code,morton_list[0].its_id);
//				printf("its tree %d\n", mpi_id);
				if(mpi_id>slave_num+1)
					return;
//				printf("its tree slave over\n");
				slave_pre_estimation(morton_list, its_queue_range, its_task_num, photon_bvh_tree,argv[1]);
			}
				
		}
	}
	athread_halt();
	MPI_Finalize();
	
	return 1;
}
