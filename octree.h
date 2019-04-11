#include "stdio.h"
#include <athread.h>
extern SLAVE_FUN(PSRS_sort_part1)(void*);
extern SLAVE_FUN(PSRS_sort_part2)(void*);
struct Intersection
{
	float p[3];
	float n[3];
};
struct Its_Queue
{
	float bound_box[6];
	int range[2];
};
struct Photon_Tree_Node
{
	float bound_box[6];
	int leaf_child;
	int children[2];
	int range[2];
};
struct BVH_tree
{
	float box[6];
	int children[2];
	int file_id;
	int range[2];
};
typedef struct
{
	struct Morton_node *M;
	int *M_S;
	int *T;
	int *S_D;
	int *S_N;
} Info_sort;
void aabbcpy(float *a, float *b){
	int i = 0;
	for (i; i < 6; i++){
		a[i] = b[i];
	}
}

int its_tree_size = 0;
int photon_tree_size = 0;
void set_photon_tree_node(struct Photon_Tree_Node *photon_tree, float *aabb, int *children, int st, int end){
	aabbcpy(photon_tree->bound_box, aabb);
	photon_tree->children[0] = children[0];
	photon_tree->children[1] = children[1];
	photon_tree->leaf_child = children[2];
	photon_tree->range[0] = st;
	photon_tree->range[1] = end;
}

int partition_i(struct Intersection *its_list, int *its_id, int st, int end, int axis)
{
	float x = its_list[st].p[axis];
	int i = st;
	int j = end + 1;
	do
	{
		do
		{
			j--;
		} while (x < its_list[j].p[axis] && j>=st);

		do
		{
			i++;
		} while (x >its_list[i].p[axis] && i<end);

		if (i < j)
		{
			struct Intersection temp = its_list[i];
			int temp_id = its_id[i];
			its_list[i] = its_list[j];
			its_id[i] = its_id[j];
			its_list[j] = temp;
			its_id[j] = temp_id;
		}
	} while (i < j);
//	swap(st, j, its_list);
	struct Intersection temp = its_list[st];
	int temp_id = its_id[st];
	its_list[st] = its_list[j];
	its_id[st] = its_id[j];
	its_list[j] = temp;
	its_id[j] = temp_id;
	return j;           // returns middle subscript  
}

void nth_element_i(struct Intersection *its_list, int *its_id, int top, int split, int bottom, int axis){
	int mid = partition_i(its_list, its_id, top, bottom, axis);
	while (mid != split)
	{
		if (split < mid){
			bottom = mid - 1;
			mid = partition_i(its_list,its_id, top, bottom, axis);
					
		}
		else{
			top = mid + 1;
			mid = partition_i(its_list,its_id,top, bottom, axis);
		}
	}
}

void get_bound_box_i(struct Intersection *list, int st, int end, float *aabb){
	int i, k;
	for (i = st; i < end; i++){
		for (k = 0; k < 3; k++){
			if (list[i].p[k] < aabb[k]|| i == st){
				aabb[k] = list[i].p[k];
			}
			if (list[i].p[k] > aabb[k + 3]|| i == st){
				aabb[k + 3] = list[i].p[k];
			}
		}
	}
}

int partition_i2(struct Intersection *its_list, int st, int end, int axis)
{
	float x = its_list[st].p[axis];
	int i = st;
	int j = end + 1;
	do
	{
		do
		{
			j--;
		} while (x < its_list[j].p[axis] && j>=st);

		do
		{
			i++;
		} while (x >its_list[i].p[axis] && i<end);

		if (i < j)
		{
			struct Intersection temp = its_list[i];
			its_list[i] = its_list[j];
			its_list[j] = temp;
		}
	} while (i < j);
//	swap(st, j, its_list);
	struct Intersection temp = its_list[st];
	its_list[st] = its_list[j];
	its_list[j] = temp;
	return j;           // returns middle subscript  
}

void nth_element_i2(struct Intersection *its_list, int top, int split, int bottom, int axis){
	int mid = partition_i2(its_list, top, bottom, axis);
	while (mid != split)
	{
		if (split < mid){
			bottom = mid - 1;
			mid = partition_i2(its_list, top, bottom, axis);
					
		}
		else{
			top = mid + 1;
			mid = partition_i2(its_list,top, bottom, axis);
		}
	}
}
int d_ab_n = 0;
void divid_its_task(struct Intersection *s_i, float aabb[6],
				float *aabb_queue, int st, int end, int mul_dep){
	if (mul_dep==0){	
		int ab_p = d_ab_n*6;
		d_ab_n++;
		int  k = 0;
		for(k = 0;k<6;k++){
			aabb_queue[ab_p + k] = aabb[k];
		}
		return;
	}
	int axis = getLargestAxis_old(&aabb[3], aabb);
	int middle = st + (end - st) / 2;
	nth_element_i2(s_i, st, middle, end, axis);
	
	float aabb_child[6] = {0};
	get_bound_box_i(s_i, st, middle, aabb_child);  //重新计算包围盒 	
	divid_its_task(s_i, aabb_child, aabb_queue, st, middle, mul_dep-1); 
	
	get_bound_box_i(s_i, middle, end, aabb_child);
	divid_its_task(s_i, aabb_child, aabb_queue, middle, end, mul_dep-1); 
	return;
}
void get_bound_box_p(struct Photon *photon_list, int st, int end, float *aabb){
	int i, k;
	for(k = 0;k<3;k++){
		aabb[k] = photon_list[st].pos[k];
		aabb[k+3] = photon_list[st].pos[k];
	}	
	for (i = st+1; i < end; i++){
		for (k = 0; k < 3; k++){		
			if (photon_list[i].pos[k] < aabb[k]){
				aabb[k] = photon_list[i].pos[k];
			}
			if (photon_list[i].pos[k] > aabb[k + 3]){
				aabb[k + 3] = photon_list[i].pos[k];
			}
		}
	}
}

int box_intsect(float *aabb1,float *aabb2){   //1包围2 返回1 否则返回0
	int i;
	for (i = 0; i < 3; i++){
		if (aabb1[i]>aabb2[i] || aabb1[3 + i] < aabb2[3 + i])
			return 0;
	}
	return 1;
}
int box_intsect2(float *aabb1, float *aabb2){   //1包围2 返回1 否则返回0
	int i;
	for (i = 0; i < 3; i++){
		if (aabb1[i]>aabb2[i] || aabb1[3 + i] < aabb2[3 + i])
			return 0;
	}
	return 1;
}
int have_intsect(float *aabb1, float *aabb2){//相交则返回1 否则返回0
	int i;
	for (i = 0; i < 3; i++){
		if (aabb1[i]>aabb2[3 + i] || aabb1[3 + i] < aabb2[i])
			return 0;
	}
	return 1;
}
int pos_aabb_intsect(float *pos, float *aabb){   //1包围2 返回1 否则返回0
	int i;
	for (i = 0; i < 3; i++){
		if (pos[i] < aabb[i] || pos[i] > aabb[3 + i])
			return 0;
	}
	return 1;
}

//广度优先遍历求出节点
float* box_merge(float *aabb1, float *aabb2){   // 返回两个盒子的 总包围盒
	int i;
	float *aabb3 = (float*)malloc(6*sizeof(float));
	for (i = 0; i < 3; i++){
		aabb3[i] = aabb1[i]<aabb2[i]?aabb1[i]:aabb2[i];
		aabb3[i+3] = aabb1[i+3]>aabb2[i+3]?aabb1[i+3]:aabb2[i+3];
	}
	return aabb3;
}
void multi_box_merge(float *aabb_list, float *res, int n){   // 返回两个盒子的 总包围盒
	int i,k;
	for(i = 0;i<n;i++){
		for (k = 0; k < 3; k++){
			if(i==0||aabb_list[i*6+k] < res[k]){
				res[k] = aabb_list[i*6+k];
			}
			if(i==0||aabb_list[i*6+3+k] > res[k+3]){
				res[k+3] = aabb_list[i*6+3+k];
			}
		}
	}
}
void set_bvh_node(struct BVH_tree *bvh_node, float *aabb, 
	int child1,int chld2, int file_id, int range1, int range2){
	aabbcpy(bvh_node->box, aabb);
	bvh_node->children[0] = child1;
	bvh_node->children[1] = chld2;
	bvh_node->file_id = file_id;
	bvh_node->range[0] = range1;
	bvh_node->range[1] = range2;
}
void get_top_tree(struct Photon *sub_tree, struct BVH_tree *bvh_tree,
	int tree_lv, float *aabb, int mpi_id, int tree_size, int photon_size){
		
	int sub_tree_index[5000];	   //存的是sub tree 的index 
	sub_tree_index[0] = 0;
	int head = 0, rear = 1;
	int base_index = mpi_id*tree_size;
	
	set_bvh_node(&bvh_tree[0], aabb, 1, 2, mpi_id, 0, photon_size);

	while (rear < tree_size){
		int cur = sub_tree_index[head];
		struct BVH_tree *bvh_node = &bvh_tree[head];
		struct Photon *node = &sub_tree[cur];
		if (isLeaf(node->flag) || node->right == 0){
			//如果sub tree 的节点是叶子节点 或者只有一个节点  则 bvh_node是叶子节点 
			bvh_node->children[0] = -1;
			bvh_node->children[1] = -1; 
			//不能直接跳过 否则乱套了  或者其实不会乱套
			bvh_tree[rear] = *bvh_node;
			sub_tree_index[rear++] = cur;
			bvh_tree[rear] = *bvh_node;
			sub_tree_index[rear++] = cur;
		}
		else{
			float aabb_temp[6];
			get_bound_box_p(sub_tree, bvh_tree[head].range[0] + 1, node->right, aabb_temp);	

			set_bvh_node(&bvh_tree[rear], aabb_temp,
				-1, -1, mpi_id, bvh_tree[head].range[0] + 1, node->right);
			bvh_node->children[0] = rear+base_index;
			sub_tree_index[rear++] = cur + 1;
			
		
			get_bound_box_p(sub_tree, node->right, bvh_tree[head].range[1], aabb_temp);	
			set_bvh_node(&bvh_tree[rear], aabb_temp,
				-1, -1, mpi_id, node->right, bvh_tree[head].range[1]);
			bvh_node->children[1] = rear+base_index;
			sub_tree_index[rear++] = node->right;
			
			
/*			float aabb_temp[6];
			aabbcpy(aabb_temp, bvh_node->box);
			int axis = node->flag;
			float midpoint = 0.5f * (aabb_temp[axis + 3] + aabb_temp[axis]);
			float temp = aabb_temp[axis + 3];
			aabb_temp[axis + 3] = midpoint;

			set_bvh_node(&bvh_tree[rear], aabb_temp,
				-1, -1, mpi_id, bvh_tree[head].range[0] + 1, node->right);
			bvh_node->children[0] = rear+base_index;
			sub_tree_index[rear++] = cur + 1;
			
		
			aabb_temp[axis + 3] = temp;
			aabb_temp[axis] = midpoint;
			set_bvh_node(&bvh_tree[rear], aabb_temp,
				-1, -1, mpi_id, node->right, bvh_tree[head].range[1]);
			bvh_node->children[1] = rear+base_index;
			sub_tree_index[rear++] = node->right;*/
			
		
		}
		head++;
	}
	int i;
	for (i = head; i < rear; i++){
		bvh_tree[i].children[0] = -1;
		bvh_tree[i].children[1] = -1;
	}
}



void construct_photon_bvh(int *subListSizes, 
			struct BVH_tree *photon_bvh_tree, int numprocs, int slave_top_tree_size){
	int i;
	int origin_size = slave_top_tree_size*numprocs;
	int st = 0, ed = origin_size;
	while(st + 1 < ed){
		float *aabb_ed; 
		if(st<origin_size){
			aabb_ed = box_merge(photon_bvh_tree[st].box,photon_bvh_tree[st+slave_top_tree_size].box);
			set_bvh_node(&photon_bvh_tree[ed],aabb_ed, 
					st, st + slave_top_tree_size, -1, -1, -1);
			st+=slave_top_tree_size*2; ed++;
		}
		else{
			aabb_ed = box_merge(photon_bvh_tree[st].box,photon_bvh_tree[st+1].box);
			set_bvh_node(&photon_bvh_tree[ed],aabb_ed, 
					st, st + 1, -1, -1, -1);	
			st+=2; ed++;
		}
//		printf("%d %d\n",st,ed);
	}
//	printf("%d \n",st);
}
int traverse_BVH_tree(struct BVH_tree *photon_bvh_tree, float *aabb, 
	int *cut,float *cut_aabb, int radius, int tree_size){
	
//	int tree_size = (numprocs-1)*8-1;
	int cut_num = 0;	
	int queue[1000];
	queue[0] = tree_size - 1;//root_index;
	struct BVH_tree *bvh_node;
	int num_current = 1;
//	printf("%d \n",tree_size);
	while (num_current != 0){
		int node_index = queue[--num_current];
		bvh_node = &photon_bvh_tree[node_index];
		
		int i;
		
//		for(i=0;i<6;i++){
//			printf("%f ",bvh_node->box[i]);
//		}
		
		if (have_intsect(aabb, bvh_node->box)){
			if(bvh_node->file_id ==-1){
				queue[num_current++] = bvh_node->children[0];
				queue[num_current++] = bvh_node->children[1];
			}
			else{//numprocs>64||
				if(bvh_node->children[0]==-1||box_intsect(aabb, bvh_node->box)){
//					printf("%d    ",node_index);
					cut[cut_num*3] = bvh_node->file_id;
					cut[cut_num*3 + 1] = bvh_node->range[0];
					cut[cut_num*3 + 2] = bvh_node->range[1];
					int i;
					for(i = 0;i<3;i++){
						cut_aabb[cut_num*6+i] = bvh_node->box[i] - radius;
						cut_aabb[cut_num*6+i+3] = bvh_node->box[i+3] + radius;
					}
					cut_num++;
				}
				else{
					queue[num_current++] = bvh_node->children[0];
					queue[num_current++] = bvh_node->children[1];
				}
			}
		}
	}
	return cut_num;
}


struct Morton_node
{
	int code;
	int its_id; //st ed
};
struct M_Tree_node
{
	int range[2]; //st ed  节点范围
	int leaf_child;   //是否为叶子  其实如果range = children就是叶子
	int children[2];  //树的子节点
	float aabb[6];   //包围盒
};
#define DEP 10
int mortonid_ph(struct Photon d, float *min, float range[3]){
	int id = 0, pos[3];
	int mask = 0;
	int i;
	pos[0] = 1000 * (d.pos[0]-min[0]) / range[0];
	pos[1] = 1000 * (d.pos[1]-min[1]) / range[1];
	pos[2] = 1000 * (d.pos[2]-min[2]) / range[2];
	for (i = 0; i<DEP + 1; i++)
	{
		mask = 1 << i;
		id += ((!!(pos[2] & mask)) * 4 + (!!(pos[1] & mask)) * 2 + (!!(pos[0] & mask)) * 1)*((int)(pow((double)8, (double)i)));
	}
	return id;
}
int mortonid_its(struct Intersection d, float *min, float range[3]){
	int id = 0, pos[3];
	int mask = 0;
	int i;
	pos[0] = 1023 * (d.p[0]-min[0]) / range[0];
	pos[1] = 1023 * (d.p[1]-min[1]) / range[1];
	pos[2] = 1023 * (d.p[2]-min[2]) / range[2];
	for (i = 0; i<DEP + 1; i++)
	{
		mask = 1 << i;
		id += ((!!(pos[2] & mask)) * 4 + (!!(pos[1] & mask)) * 2 + (!!(pos[0] & mask)) * 1)*((int)(pow((double)8, (double)i)));
	}
	return id;
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


int PSRS_sort(struct Morton_node *morton_list, struct Morton_node *new_morton_list, int size){
	int j,i;

	int sort_task_num = size / 64;
	int sample_num = sort_task_num / 64;
	struct Morton_node *main_element;// = (struct Morton_node *)malloc(5000 * sizeof(struct Morton_node));
	int main_num = 0;
	
	athread_spawn(PSRS_sort_part1, morton_list);
	athread_join();
	
	main_element = &morton_list[size];
	quickSort(main_element, 0, 64*64-1);   //代表元素排序
	int main_element_sample[65];  
	int thread_tag[65*65]; 
	for (i = 0; i < 64; i++){
		main_element_sample[i] = main_element[i * 64].code;    //代表元素正则采样
	}
	main_element_sample[64] = 2147483647;
	int *sorted_data = (int *)malloc(size * 2 * sizeof(int));
	int thread_data_num[64];
	
	for(i = 0;i<65*65;i++){
		thread_tag[i] = 0;
	}
	
/*	for(i = 0;i<65;i++){
		thread_tag[i*65] = i*sort_task_num;
		int current = 0;
		for (j = 1; j < 65; j++){
			int current_main_tag = main_element_sample[j];    //小于这个tag的
			while (morton_list[i*sort_task_num+current].code < current_main_tag && current < sort_task_num){
				current++;
			}
			thread_tag[i*65 + j] = current + i*sort_task_num;
		}
	}*/
	
	Info_sort info;
	info.M = morton_list;
	info.M_S = &main_element_sample[0];
	info.T = &thread_tag[0];
	info.S_D = &sorted_data[0];
	info.S_N = &thread_data_num[0];
	athread_spawn(PSRS_sort_part2, &info);
	athread_join();
	
	int date_range = sort_task_num*2;
	for (i = 0; i < 64; i++){
		thread_data_num[i] = 0;
		int thread_current[64];
		int dead = 0;
		for (j = 0; j < 64; j++){
			//第i个线程 取所有第i个段
			if (thread_tag[j*65 + i + 1] - thread_tag[j*65 + i] > 0){
				thread_current[j] = thread_tag[j*65 + i];   
			}
			else{
				thread_current[j] = -1;
				dead++;
			}
		}
		while (dead < 64){
			int min = 2147483647, min_index = 0, min_j = 0;
			
			//每次从里边找个最小的放进去
			for (j = 0; j < 64; j++){
				if (thread_current[j] == -1)
					continue;
				if (morton_list[thread_current[j]].code < min){
					min = morton_list[thread_current[j]].code;
					min_index = thread_current[j];
					min_j = j;
				}
			}
			sorted_data[i * date_range + thread_data_num[i]++] = min_index;
			if (thread_current[min_j] + 1 < thread_tag[min_j*65 + i + 1]) //开区间
				thread_current[min_j]++;
			else{
				thread_current[min_j] = -1;
				dead++;
			}
		}
	}
	int all_num = 0;
	//排序正确?
	for (i = 0; i < 64; i++){
		//thread_data[all_num]
		for (j = 0; j < thread_data_num[i]; j++){
			sorted_data[all_num] = sorted_data[i * date_range + j];
			new_morton_list[all_num++] = morton_list[sorted_data[i * date_range + j]];
			if (new_morton_list[all_num - 1].code<new_morton_list[all_num - 2].code){
				printf("sort  error%d \n", all_num);
			}
		}
		int aaa = 0;
	}
	free(main_element);
	free(sorted_data);
	return all_num;
}

void get_bound_box_morton(struct Morton_node *morton_list, 
				struct Intersection *list, int num, float *aabb){
	int i, k;
	for (k = 0; k < 3; k++){					
		aabb[k] = 0;
		aabb[k + 3] = 0;	
	}
	for (i = 0; i < num; i++){
		struct Intersection *its_node = &list[morton_list[i].its_id];
//		printf("%d  %d  %d   ", num, i, morton_list[i].its_id);
//		printf("%f  %f   %f\n",its_node->p[0],its_node->p[1],its_node->p[2]);
		
		for (k = 0; k < 3; k++){			
			if (its_node->p[k] < aabb[k]|| i == 0){
				aabb[k] = its_node->p[k];
			}
			if (its_node->p[k] > aabb[k + 3]|| i == 0){
				aabb[k + 3] = its_node->p[k];
			}
		}
	}
	for(i = 0; i<3; i++){
		aabb[i] -= 0.000002;
		aabb[i+3] += 0.000002;
	}
}
void construct_morton_queue(struct Intersection *its_list, struct Morton_node *morton_list, 
						int min, int max,float *queue_aabb, int *queue_range,
						int its_size, int *num, int *head, int *rear, int *change_point){
	
	if (morton_list[(*head)].code > max||morton_list[(*rear)].code < min||(*head)==(*rear))
		return;
	if (morton_list[(*rear)].code >= max){  //code 小 可以组成一个节点
		
		int n = 0, st = (*head);
//		printf("code0 %d  %d %d",min,max, (*num));		
//		printf("code %d  %d %d",morton_list[(*head)].code, morton_list[(*rear)].code, max);
		while(morton_list[(*head)].code <= max&&(*head)<its_size){
			n++;
			(*head)++;// = (?its_size-1:(*head)+1;
			(*rear) = ((*rear)==its_size)?its_size:(*rear)+1;
//			printf("code2 %d %d %d  %d %d", (*head),(*rear),morton_list[(*head)].code, morton_list[(*rear)].code, max);
		}
		
		int wirte_pos = 6*(*num);
		get_bound_box_morton(&morton_list[st], its_list, n, &queue_aabb[wirte_pos]);		
		queue_range[(*num)] = n;
		
		if((*change_point) != 0&&(*num)>(*change_point)){
			(*rear) = (*head) + ((*rear) - (*head))/2;    //Reduced 
			(*change_point) = 0;
		}
			
//		int j;		
//		for(j = 0;j<6;j++)
//			printf("%f ",queue_aabb[wirte_pos + j]);
//		printf("nn %d\n",n);
		
		(*num)++;
		return;

	}
	else{
		construct_morton_queue(its_list, morton_list, min, min + (max - min) / 2, 
								queue_aabb, queue_range, its_size, num, head, rear, change_point);
							
		construct_morton_queue(its_list, morton_list, min + (max - min) / 2, max, 
								queue_aabb, queue_range, its_size, num, head, rear, change_point);					
	}
	return;
}

