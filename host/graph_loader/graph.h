#ifndef __GRAPH_H__
#define __GRAPH_H__

#include "host_config.h"
#include <vector>
#include <fstream>
#include <iostream>
#include <iterator>
#include <cstdlib>
#include <sstream>
#include <cmath>
#include <map>
#include <algorithm>

#define HERE do {std::cout << "File: " << __FILE__ << " Line: " << __LINE__ << std::endl;} while(0)

class Vertex {
    public:
        int idx;
        int inDeg;
        int outDeg;

        std::vector<int> inVid;
        std::vector<int> outVid;

        explicit Vertex(int _idx) {
            idx = _idx;
        }

        ~Vertex(){
            // Nothing is done here.
        }

};

class Graph{
    public:
        int vertexNum;
        int edgeNum;
        std::vector<Vertex*> vertices; 

        Graph(const std::string &fName);
        ~Graph(){
			for(int i = 0; i < vertexNum; i++){
				delete vertices[i];
			}
		};
        void getRandomStartIndices(std::vector<int> &startIndices);
        void getStat();

    private:
        bool isUgraph;
        int getMaxIdx(const std::vector<std::vector<int>> &data);
        int getMinIdx(const std::vector<std::vector<int>> &data);
        void loadFile(
				const std::string& fName,
                std::vector<std::vector<int>> &data
				);

};


/*
图为 0→2, 0→3, 1→2
顶点ID:  0      1      2
出边:    2,3    2      (没有)
入边:    (没有) (没有)  0,1

出边结构（Outgoing Edges）
rpao = [0, 2, 3, 3]  // 前缀和，表示累计的出边数
ciao = [2, 3, 2]     // 存储所有出边的目标顶点

入边结构（Incoming Edges）  
rpai = [0, 0, 0, 2]  // 前缀和，表示累计的入边数
ciai = [0, 1]        // 存储所有入边的源顶点
*/
// 用于储存图结构
class CSR{
    public:
        // 顶点总数
        const int vertexNum;
        // 边总数
        const int edgeNum;
        // rpao[i] 是前 i 个顶点的出边总数
        std::vector<int> rpao; 
        // 存储所有顶点的出边目标
        std::vector<int> ciao;
        // rpai[i] 是前 i 个顶点的入边总数
        std::vector<int> rpai;
        // 存储所有顶点的入边源
        std::vector<int> ciai;
        // 边属性
        std::vector<prop_t> eProps;
        // 顶点属性
        std::vector<prop_t> vProps;

        // The CSR is constructed based on the simple graph
        explicit CSR(const Graph &g);
        int save2File(const std::string & fName);
		~CSR();
};

class CSR_BLOCK{
	public:
		const int cordx;
		const int cordy;
		int vertexNum;
		int edgeNum;
		int srcStart;
		int srcEnd;
		int sinkStart;
		int sinkEnd;
		std::vector<int> rpa;
		std::vector<int> cia;
		std::vector<prop_t> eProps;
		explicit CSR_BLOCK(const int _cordx, const int _cordy, CSR* csr);
		~CSR_BLOCK(){};
};

#endif
