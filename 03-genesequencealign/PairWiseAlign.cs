using System;
using System.Collections.Generic;
using System.Text;

namespace GeneticsLab
{
    enum operation { Top, Left, Sub, Match, first};
    class Node
    {
        public Node prev { get; set; }
        public operation op { get; set; }
        public int value { get; set; }
        public Node(Node Prev, operation op, int value)
        {
            this.prev = prev;
            this.op = op;
            this.value = value;
        }
    }
    class PairWiseAlign
    {
        int MaxCharactersToAlign;

        public PairWiseAlign()
        {
            // Default is to align only 5000 characters in each sequence.
            this.MaxCharactersToAlign = 5000;
        }

        public PairWiseAlign(int len)
        {
            // Alternatively, we can use an different length; typically used with the banded option checked.
            this.MaxCharactersToAlign = len;
        }

        /// <summary>
        /// this is the function you implement.
        /// </summary>
        /// <param name="sequenceA">the first sequence</param>
        /// <param name="sequenceB">the second sequence, may have length not equal to the length of the first seq.</param>
        /// <param name="banded">true if alignment should be band limited.</param>
        /// <returns>the alignment score and the alignment (in a Result object) for sequenceA and sequenceB.  The calling function places the result in the dispay appropriately.
        /// 
        public ResultTable.Result Align_And_Extract(GeneSequence sequenceA, GeneSequence sequenceB, bool banded)
        {
            ResultTable.Result result = new ResultTable.Result();
            int score;                                                       // place your computed alignment score here
            string[] alignment = new string[2];                              // place your two computed alignments here

            int sizeA = MaxCharactersToAlign;
            int sizeB = MaxCharactersToAlign;

            if (sequenceA.Sequence.Length < MaxCharactersToAlign) sizeA = sequenceA.Sequence.Length;
            if (sequenceB.Sequence.Length < MaxCharactersToAlign) sizeB = sequenceB.Sequence.Length;

            if (sequenceA == sequenceB)
            {
                score = -3 *sizeA;
                alignment[0] = sequenceA.Sequence;
                alignment[1] = sequenceB.Sequence;
                // ***************************************************************************************


                result.Update(score, alignment[0], alignment[1]);                  // bundling your results into the right object type 
                return (result);
            }

            int[,] nodes = new int[sizeA+1,sizeB+1];
            operation[,] prevs = new operation[sizeA + 1, sizeB + 1];

          
            for (int i = 0; i < sizeA + 1; i++)
            {
                for (int j = 0; j < sizeB + 1; j++)
                {
                    if (banded && j > i + 3) break;
                    if (banded && j < i - 3) continue;

                    if (i == 0)
                    {
                        if (j == 0)
                        {

                            nodes[i, j] = 0;
                        }
                        else
                        {
                            nodes[i, j] = nodes[i, j - 1] + 5;
                        }
                    }
                    else if (j == 0)
                    {
                        nodes[i, j] = nodes[i - 1, j] + 5;
                    }
                    else
                    {

                        int min_match;
                        int min_InDel_top;
                        int min_InDel_left;
                        int min_sub;
                        if (sequenceA.Sequence[i - 1] == sequenceB.Sequence[j - 1])
                        {
                            min_match = nodes[i - 1, j - 1] - 3;

                        }
                        else
                        {
                            min_match = int.MaxValue;
                        }
                        if (banded && j > i - 1 + 3)
                        {
                            min_InDel_top = int.MaxValue;
                        }
                        else
                        {
                            min_InDel_top = nodes[i - 1, j] + 5;
                        }
                        if (banded && j - 1 < i - 3)
                        {
                            min_InDel_left = int.MaxValue;
                        }
                        else
                        {
                            min_InDel_left = nodes[i, j - 1] + 5;
                        }


                        min_sub = nodes[i - 1, j - 1] + 1;

                        if (min_match <= min_InDel_top && min_match <= min_InDel_left && min_match <= min_sub)
                        {
                            nodes[i, j] = min_match;
                            prevs[i, j] = operation.Match;
                        }
                        else if (min_sub < min_match && min_sub <= min_InDel_top && min_sub <= min_InDel_left)
                        {
                            nodes[i, j] = min_sub;
                            prevs[i, j] = operation.Sub;
                        }
                        else if (min_InDel_top < min_match && min_InDel_top < min_sub && min_InDel_top < min_InDel_left)
                        {
                            nodes[i, j] = min_InDel_top;
                            prevs[i, j] = operation.Top;
                        }
                        else
                        {
                            nodes[i, j] = min_InDel_left;
                            prevs[i, j] = operation.Left;
                        }

                    }
                    
                }
            }
            
            int countA = sizeA;
            int countB = sizeB;
            if(banded && Math.Abs(sizeA - sizeB) > 7)
            {
                
                alignment[0] = "No Alignment Possible";
                alignment[1] = "No Alignment Possible";
                score = int.MaxValue;
            }
            else
            {
                score = nodes[sizeA, sizeB];
                while (countA > 0 && countB > 0)
                {
                    if (prevs[countA, countB] == operation.Match || prevs[countA, countB] == operation.Sub)
                    {
                        alignment[0] = sequenceA.Sequence[countA - 1] + alignment[0];
                        alignment[1] = sequenceB.Sequence[countB - 1] + alignment[1];
                        countA--;
                        countB--;
                    }
                    else if (prevs[countA, countB] == operation.Top)
                    {
                        alignment[0] = sequenceA.Sequence[countA - 1] + alignment[0];
                        alignment[1] = '-' + alignment[1];
                        countA--;
                    }
                    else if (prevs[countA, countB] == operation.Left)
                    {
                        alignment[0] = '-' + alignment[0];
                        alignment[1] = sequenceB.Sequence[countB - 1] + alignment[1];
                        countB--;
                    }

                }
            }
            

            /*for (int i = 0; i < sizeA+1; i++)
            {
                for (int j = 0; j < sizeB+1; j++)
                {
                    Console.Write(nodes[i, j] + " "+prevs[i,j]);
                  
                }
                Console.WriteLine();
            }*/

            result.Update(score,alignment[0],alignment[1]);  // bundling your results into the right object type 
            return(result);
        }
    }
}
