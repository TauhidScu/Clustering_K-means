#include <stdlib.h>
#include <vector>
#include <string>
#include <queue>

const std::string g_FileName = "UniformData.txt";
const int g_NumCluster = 8;

struct HuffmanTreeNode
{
	int NodeIndex, Left, Right, Parent;
};

struct QueueElement
{
	int Frequency;
	int NodeIndex;

	friend bool operator<(const QueueElement& vLeft, const QueueElement& vRight) { return vLeft.Frequency > vRight.Frequency; }
};

struct CodeTableElement
{
	int Data;
	int Frequency;
	std::vector<bool> Code;
};

void readData(std::vector<std::pair<float, float>>& voData);
void executeKMeans(std::vector<std::pair<float, float>>& vInputData, std::vector<std::pair<float, float>>& voClusterCenter, std::vector<int>& voClusterIndex);
void initClusterCenter(std::vector<std::pair<float, float>>& vInputData, std::vector<std::pair<float, float>>& voClusterCenter);
void updateClusterIndex(std::vector<std::pair<float, float>>& vInputData, std::vector<std::pair<float, float>>& vClusterCenter, std::vector<int>& voClusterIndex);
void updateClusterCenter(std::vector<std::pair<float, float>>& vInputData, std::vector<int>& vClusterIndex, std::vector<std::pair<float, float>>& voClusterCenter);
void countNumDataInCluster(std::vector<int>& vClusterIndex, std::vector<int>& voNumDataInCluster);
void executeHuffmanEncoding(const std::vector<int>& vInput, std::vector<bool>& voOutput);
void generateFrequencyTable(const std::vector<int>& vInput, std::vector<CodeTableElement>& voCodeTable);
void initPriorityQueueAndTree(std::vector<CodeTableElement>& vCodeTable, std::priority_queue<QueueElement>& voPriorityQueue, std::vector<HuffmanTreeNode>& voHuffmanTree);
void generateFrequencyTable(const std::vector<int>& vInput, std::vector<CodeTableElement>& voCodeTable);
void generateHuffmanTree(std::priority_queue<QueueElement>& vPriorityQueue, std::vector<HuffmanTreeNode>& voHuffmanTree);
void generateHuffmanCode(std::vector<HuffmanTreeNode>& vHuffmanTree, std::vector<CodeTableElement>& voCodeTable);
void generateEncodedBitStream(const std::vector<int>& vInputData, std::vector<CodeTableElement>& vCodeTable, std::vector<bool>& voOutput);

float computeError(std::vector<std::pair<float, float>>& vInputData, std::vector<int>& vClusterIndex, std::vector<std::pair<float, float>>& vClusterCenter);
float computeDistance(std::pair<float, float>& v1, std::pair<float, float>& v2);

int assignClusterIndex(std::pair<float, float>& vInput, std::vector<std::pair<float, float>> vClusterCenter);

void initPriorityQueueAndTree(std::vector<CodeTableElement>& vCodeTable, std::priority_queue<QueueElement>& voPriorityQueue, std::vector<HuffmanTreeNode>& voHuffmanTree)
{
	QueueElement q;
	HuffmanTreeNode t;

	for (int i=0; i<vCodeTable.size(); i++)
	{
		q.Frequency = vCodeTable[i].Frequency;
		q.NodeIndex = i;
		voPriorityQueue.push(q);

		t.Left = t.Right = t.Parent = -1;
		t.NodeIndex = i;
		voHuffmanTree.push_back(t);
	}
}

void generateFrequencyTable(const std::vector<int>& vInput, std::vector<CodeTableElement>& voCodeTable)
{
	for(int i = 0; i < vInput.size(); i++)
	{	
		int j;
		for(j = 0; j < voCodeTable.size(); j++)
		{		
			if(voCodeTable[j].Data == vInput[i])
			{
				voCodeTable[j].Frequency++;			
				break;
			}
		}
		if (j == voCodeTable.size())
		{
			CodeTableElement t;
			t.Data = vInput[i];
			t.Frequency = 1;
			voCodeTable.push_back(t);
		}
	} 

	_ASSERT(voCodeTable.size() == g_NumCluster);
	int Total = 0;
	for (int i=0; i<voCodeTable.size(); i++) Total += voCodeTable[i].Frequency;
	_ASSERT(Total == vInput.size());
}

void generateHuffmanTree(std::priority_queue<QueueElement>& vPriorityQueue, std::vector<HuffmanTreeNode>& voHuffmanTree)
{
	QueueElement Left, Right, t;
	HuffmanTreeNode Node;

	while (vPriorityQueue.size() > 1)
	{
        Left=vPriorityQueue.top();
		vPriorityQueue.pop();
		Right=vPriorityQueue.top();
		vPriorityQueue.pop();

		t.Frequency= Left.Frequency+Right.Frequency;
		t.NodeIndex=voHuffmanTree.size();
		vPriorityQueue.push(t);

		Node.Left = Left.NodeIndex;
		Node.Right = Right.NodeIndex;
		Node.Parent = -1;
		Node.NodeIndex = t.NodeIndex;
	
		voHuffmanTree[Left.NodeIndex].Parent = Node.NodeIndex;
		voHuffmanTree[Right.NodeIndex].Parent = Node.NodeIndex;

		voHuffmanTree.push_back(Node);
	}
}

void generateHuffmanCode(std::vector<HuffmanTreeNode>& vHuffmanTree, std::vector<CodeTableElement>& voCodeTable)
{
	for (int i = 0; i < 8; i++)
	{
		std::vector<bool> ReverseCode;

		HuffmanTreeNode Node = vHuffmanTree[i];
		HuffmanTreeNode ParentNode;

		while (Node.Parent != -1)
		{
			ParentNode = vHuffmanTree[Node.Parent];
			ReverseCode.push_back(ParentNode.Left == Node.NodeIndex);
			Node = ParentNode;
		}

		for (int j = ReverseCode.size() - 1; j >= 0; j--)
		{
			voCodeTable[i].Code.push_back(ReverseCode[j]);
		}
	}
}

void generateEncodedBitStream(const std::vector<int>& vInputData, std::vector<CodeTableElement>& vCodeTable, std::vector<bool>& voOutput)
{
	for (int i = 0; i < vInputData.size(); i++)
	{
		for (int j = 0; j < 8; j++)
			if (vInputData[i] == vCodeTable[j].Data)
				for (bool b : vCodeTable[j].Code)
					voOutput.push_back(b);
	}
}

void executeHuffmanEncoding(const std::vector<int>& vInputData, std::vector<bool>& voOutput)
{
	std::vector<HuffmanTreeNode>  HuffmanTree;
	std::vector<CodeTableElement> CodeTable;

	generateFrequencyTable(vInputData, CodeTable);

	std::priority_queue<QueueElement> PriorityQueue;

	initPriorityQueueAndTree(CodeTable, PriorityQueue, HuffmanTree);
	generateHuffmanTree(PriorityQueue, HuffmanTree);

	generateHuffmanCode(HuffmanTree, CodeTable);
	generateEncodedBitStream(vInputData, CodeTable, voOutput);
}

int assignClusterIndex(std::pair<float, float>& vInput, std::vector<std::pair<float, float>> vClusterCenter)
{
	float min_distance = FLT_MAX;
	int result;
	for(int i=0;i<8;i++){

		float distance = computeDistance (vInput,vClusterCenter[i]);
		if(distance<min_distance){
			min_distance = distance;
			result = i;
		}
	}
	return result;
}

float computeDistance(std::pair<float, float>& v1, std::pair<float, float>& v2)
{
	return sqrt((v1.first - v2.first) * (v1.first - v2.first)
		+(v1.second - v2.second) * (v1.second - v2.second));
}

void updateClusterIndex(std::vector<std::pair<float, float>>& vInputData, std::vector<std::pair<float, float>>& vClusterCenter, std::vector<int>& voClusterIndex)
{
	for (int i=0; i<vInputData.size(); i++)
	{
		voClusterIndex[i] = assignClusterIndex(vInputData[i], vClusterCenter);
	}
}

void countNumDataInCluster(std::vector<int>& vClusterIndex, std::vector<int>& voNumDataInCluster)
{
	for(int i = 0; i<vClusterIndex.size(); i++)
	{
		int ClusterIndex = vClusterIndex[i];
		voNumDataInCluster[ClusterIndex]++;
	}
}

void updateClusterCenter(std::vector<std::pair<float, float>>& vInputData, std::vector<int>& vClusterIndex, std::vector<std::pair<float, float>>& voClusterCenter)
{
	std::vector<int> NumDataInCluster;
	std::vector<std::pair<float, float>> NewClusterCenter;
	
	for (int i=0; i<g_NumCluster; i++) 
	{
		NumDataInCluster.push_back(0);
		NewClusterCenter.push_back(std::make_pair(0,0));
	}

	countNumDataInCluster(vClusterIndex, NumDataInCluster);

	for (int i=0; i<vInputData.size(); i++)
	{
		NewClusterCenter[vClusterIndex[i]].first += vInputData[i].first ;
		NewClusterCenter[vClusterIndex[i]].second += vInputData[i].second;
	}

	for (int i=0; i<g_NumCluster; i++)
	{
		NewClusterCenter[i].first /= NumDataInCluster[i];
		NewClusterCenter[i].second /= NumDataInCluster[i];
	}
}

float computeError(std::vector<std::pair<float, float>>& vInputData, std::vector<int>& vClusterIndex, std::vector<std::pair<float, float>>& vClusterCenter)
{
	float Error = 0;

	for (int i=0; i<vInputData.size(); i++)
	{
		Error += computeDistance(vInputData[i], vClusterCenter[vClusterIndex[i]]);
	}

	Error /= vInputData.size();
	return Error;
}

void readData(std::vector<std::pair<float, float>>& voData)
{
	FILE *fp = fopen(g_FileName.c_str(),"r");
	int pairNum;

	fscanf(fp,"%d\n", &pairNum);

	for(int i=0;i<pairNum;i++)
	{
		float weight,height;
		fscanf(fp,"%f %f\n",&weight, &height);
		voData.push_back(std::make_pair(weight,height));
	}
}

void initClusterCenter(std::vector<std::pair<float, float>>& vInputData, std::vector<std::pair<float, float>>& voClusterCenter)
{
	voClusterCenter.resize(8);
	for(int i = 0; i < 8 ; i ++) 
		voClusterCenter[i] = vInputData[i];
}

void executeKMeans(std::vector<std::pair<float, float>>& vInputData, std::vector<std::pair<float, float>>& voClusterCenter, std::vector<int>& voClusterIndex)
{
	initClusterCenter(vInputData, voClusterCenter);

	voClusterIndex.resize(vInputData.size());

	float LastError = FLT_MAX;
	float CurrentError;
	int NumIteration = 0;

	while (NumIteration < 100)
	{
		updateClusterIndex(vInputData, voClusterCenter, voClusterIndex);
		updateClusterCenter(vInputData, voClusterIndex, voClusterCenter);

		CurrentError = computeError(vInputData, voClusterIndex, voClusterCenter);
		if (fabs(CurrentError - LastError) < 0.0001) break;
		LastError = CurrentError;
		NumIteration++;
	}
}

int main()
{
	std::vector<std::pair<float, float>> InputData;
	std::vector<std::pair<float, float>> ClusterCenter;
	std::vector<int> ClusterIndex;
	std::vector<bool> BitStream;

	readData(InputData);

	executeKMeans(InputData, ClusterCenter, ClusterIndex);
	executeHuffmanEncoding(ClusterIndex, BitStream); 
	return 0;
}
