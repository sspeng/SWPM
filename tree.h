#include "stdio.h"
enum {
	ELeafFlag = 0x10,
	EAxisMask = 0x0F
};
struct Photon
{
	float pos[3];
	unsigned int right;
	unsigned char data[8];
	unsigned short dep;
	char flag;
};
struct cPhoton
{
	float pos[3];
	unsigned int children[2];//0 start 1 num
	unsigned char data[8];
	unsigned short dep;
	char flag;
};
inline int isLeaf(char flags)  { return flags & (unsigned char)ELeafFlag; }
inline unsigned short getAxis(char flags)  { return flags & (unsigned char)EAxisMask; }
int getLargestAxis_old(float aabb_max[3], float aabb_min[3]){
	float d[3];
	d[0] = aabb_max[0] - aabb_min[0];
	d[1] = aabb_max[1] - aabb_min[1];
	d[2] = aabb_max[2] - aabb_min[2];
	int largest = 0;
	int i;
	for (i = 1; i<3; ++i)
	if (d[i] > d[largest])
		largest = i;
	return largest;
}
int getLargestAxis(float aabb_box[6]){
	float d[3];
	d[0] = aabb_box[3] - aabb_box[0];
	d[1] = aabb_box[4] - aabb_box[1];
	d[2] = aabb_box[5] - aabb_box[2];
	int largest = 0;
	int i;
	for (i = 1; i<3; ++i)
		if (d[i] > d[largest])
			largest = i;
	return largest;
}
int count_if(int rangeStart, int rangeEnd, int axis, float midpoint, struct Photon *m_nodes){
	int num = 0;
	int i;
	for (i = rangeStart; i < rangeEnd; i++){
		if (m_nodes[i].pos[axis] <= midpoint){
			num++;
		}
	}
	return num;
}
void swap(int a, int b, struct Photon *m_nodes){
	struct Photon temp = m_nodes[a];
	m_nodes[a] = m_nodes[b];
	m_nodes[b] = temp;
}
void swap_int(int *a, int *b){
	int temp = *a;
	*a = *b;
	*b = temp;
}

//////////////////
int partition_x(struct Photon *sub_tree, int top, int bottom,int axis)
{
	float x = sub_tree[top].pos[axis];
	int i = top;
	int j = bottom + 1;
	do
	{
		do
		{
			j--;
		} while (x < sub_tree[j].pos[axis] && j>=top);

		do
		{
			i++;
		} while (x >sub_tree[i].pos[axis] && i<bottom);

		if (i < j)
		{
			struct Photon temp = sub_tree[i];
			sub_tree[i] = sub_tree[j];
			sub_tree[j] = temp;
		}
	} while (i < j);
	swap(top, j, sub_tree);
	return j;           // returns middle subscript  
}

void nth_element_x2(struct Photon *sub_tree, int top, int split, int bottom, int axis){
	int mid = partition_x(sub_tree, top, bottom, axis);
	while (mid != split)
	{
		if (split < mid){
			bottom = mid - 1;
			mid = partition_x(sub_tree, top, bottom, axis);
					
		}
		else{
			top = mid + 1;
			mid = partition_x(sub_tree,top, bottom, axis);
		}
	}
}
void construct_photon_tree_xx(int dep, struct Photon *sub_tree, int start, int end, float sub_max[3], float sub_min[3]){
	if (end - start < 1) {
		/* Create a leaf node */
		sub_tree[start].flag = (unsigned char)0x10;
		return;
	}
	
	int axis = 0;
	int split;
	axis = getLargestAxis_old(sub_max, sub_min);

	float midpoint = 0.5f * (sub_max[axis] + sub_min[axis]);

	size_t nLT = count_if(start, end, axis, midpoint, sub_tree);
//	printf("%d ",dep);
	split = start + nLT;
	if(split>end||split<start)
		printf("split num error %d\n",split);
	if (split == start)
		++split;
	else if (split == end)
		--split;
	nth_element_x2(sub_tree, start, split, end, axis);

	sub_tree[split].flag = 0x00;
	sub_tree[split].flag = sub_tree[split].flag | axis;

	if (split != end)
		sub_tree[split].right = split + 1;
	else
		sub_tree[split].right = 0;

	/*splitNode.setLeftIndex((IndexType)(rangeStart - base),
	(IndexType)(rangeStart + 1 - base));*/

	swap(start, split,sub_tree);

	/* Recursively build the children */
	float temp = sub_max[axis],
		splitPos = sub_tree[start].pos[axis];
	sub_max[axis] = splitPos;
//	if (split == 32)
//		int iiii = 9;
	construct_photon_tree_xx(dep + 1, sub_tree, start + 1, split, sub_max, sub_min);
	sub_max[axis] = temp;

//	if (split == 71224 && end == 71227)
//		int yyy = 0;

	if (split != end) {
		temp = sub_min[axis];
		sub_min[axis] = splitPos;
		construct_photon_tree_xx(dep + 1, sub_tree, split + 1, end, sub_max, sub_min);
		sub_min[axis] = temp;
	}
}

//直接层次遍历得了
void photon_to_cphoton(struct Photon a, struct cPhoton *b){
	int i;
	for (i = 0; i < 3; i++){
		b->pos[i] = a.pos[i];
	}
	for (i = 0; i < 8; i++){
		b->data[i] = a.data[i];
	}
	b->dep = a.dep;
	b->flag = a.flag;
}
int new_tree_num = 0;
int construct_chunk(struct Photon *sub_tree, struct cPhoton *new_tree, int *queue, int *new_queue, int *queue_current, int node, int chunk_level){
	int wide_queue[2000];
	wide_queue[0] = node;
	int next_wide_queue[2000];
	int node_num = 1;
	int level = 0;
	int start_node = new_tree_num;
	int i;
	while (level<chunk_level&&node_num>0){
		int next_num = 0;
		int have_leaf = 0;
		int one_child = 0;
		for (i = 0; i < node_num; i++){
			int current_node = wide_queue[i];
			struct Photon *node = &sub_tree[current_node];
			
			photon_to_cphoton(*node, &new_tree[new_tree_num++]);
			if (isLeaf(node->flag)){
				have_leaf++;
				new_tree[new_tree_num - 1].children[0] = -1;
				new_tree[new_tree_num - 1].children[1] = -1;
			}else if (level != chunk_level - 1){
				new_tree[new_tree_num - 1].children[0] = new_tree_num - 1 + (i - have_leaf) * 2 + node_num - i - one_child;
				next_wide_queue[next_num++] = current_node + 1;
				if (!node->right == 0){
					new_tree[new_tree_num - 1].children[1] = new_tree_num - 1 + (i - have_leaf) * 2 + node_num - i + 1 - one_child;
					next_wide_queue[next_num++] = node->right;
					
				}
				else{
					new_tree[new_tree_num - 1].children[1] = -1;
					one_child++;
				}
				
			}
		}
		level++;
		if (level == chunk_level)
			break;
		for (i = 0; i < next_num; i++){
			wide_queue[i] = next_wide_queue[i];
		}
		node_num = next_num;
	}
	for (i = 0; i < node_num; i++){
		new_queue[(*queue_current)] = new_tree_num - node_num + i;
		queue[(*queue_current)++] = wide_queue[i];
	}
	return start_node;
}

void construct_photon_tree_chunk(int dep, struct Photon *sub_tree,struct cPhoton *new_tree, int chunk_level, int branch){
	int queue[200000];
	int new_queue[200000];
	new_queue[0] = 0;
	queue[0] = 0;
	int current = 0;
	new_tree_num = 0;
	construct_chunk(sub_tree, new_tree, queue, new_queue, &current, 0, chunk_level);
	while (current > 0){
		int node_num = queue[--current];
		int new_node_num = new_queue[current];
		struct Photon *node = &sub_tree[node_num];
		struct cPhoton *new_node = &new_tree[new_node_num];
		if (!isLeaf(node->flag)) {
			int children[2];
			new_node->children[0] = construct_chunk(sub_tree, new_tree, queue, new_queue, &current, node_num + 1, chunk_level);
			if (!node->right == 0)
				new_node->children[1] = construct_chunk(sub_tree, new_tree, queue, new_queue, &current, node->right, chunk_level);
			else
				new_node->children[1] = -1;
		}
		else{
			new_node->children[0] = -1;
			new_node->children[1] = -1;
		}
	}
}




