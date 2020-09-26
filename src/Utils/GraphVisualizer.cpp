//
// Created by sylwester on 5/13/19.
//

#include <DataStructures/FAU.h>
#include <Utils/GraphVisualizer.h>
#include <Utils/MyUtils.h>
#include <Global.h>
#include <cassert>

#include "Utils/GraphVisualizer.h"


void GraphVisualizer::writeInGraphvizFormat( Graph * G, vector<Contig*> & contigs ) {

    ADD_REVCOMP_CONTIGS = false;

    if(ADD_REVCOMP_CONTIGS){
        int S = contigs.size();
        contigs.resize(S);
        cerr << "smallest contig has size " << contigs.back()->size() << endl;
        for(int i=0; i<S; i++) contigs.push_back( MyUtils::getCompRevContig( contigs[i], &Global::READS ) );
    }



    cerr << "Writing graph from contigs in graphviz format" << endl;

    out.open( Params::outStreamFileName + "_mopp" + to_string(Params::MAX_OFFSET_PARALLEL_PATHS)
                  + "_modb" + to_string(Params::MAX_OFFSET_DANGLING_BRANCHES) + "_rsoe" + to_string(Params::REMOVE_SMALL_OVERLAP_EDGES_MIN_OVERLAP) + "-" +
                  to_string(Params::REMOVE_SMALL_OVERLAP_EDGES_NUMBER_TO_RETAIN) + ".gv");

    out << "digraph G{" << endl;

        createRelevantNodes( G, contigs );

        visualizeGraph( G, contigs );

        visualizeContigs(contigs);
    out << "}" << endl;
    out.close();

    cerr << "writing in graphviz format finished" << endl;

   if(ADD_REVCOMP_CONTIGS) contigs.resize( contigs.size() >> 1 );
}

void GraphVisualizer::visualizeContigs(vector<Contig *> &contigs) {
    cerr << "Craeting contig graph" << endl;
    int counter = 0;

    int B = contigs.size();
    if( ADD_REVCOMP_CONTIGS ) B >>= 1;

    for( int i=0; i<B; i++ ){
        visualizeContig(contigs[i], colors[ counter++ % C ]  );
    }

    for( int i=B; i<contigs.size(); i++ ){ // here i visualize revcomp contigs
        visualizeContig(contigs[i], "black"  );
    }

    cerr << "Contig graph created" << endl << endl;
}

void GraphVisualizer::visualizeContig(Contig *ctg, string color) {
    auto path = ctg->getContainedReads();

    if( color == "black" ){ // here i mark a revcomp contig
        out << "\t" << path[0].first->getId() << " -> " << path.back().first->getId() <<
            " [label=" << ctg->size() - path.back().first->size() << ", color=" << color << ", penwidth=" << 5 << "];" << endl;
        return;
    }

    int a = 0;
    int b = 1;
    int offset = 0;

    while( b < path.size() ){
        while( relevantNodes.count( path[b].first->getId() ) == 0 ){
            offset += path[b].second;
            b++;
        }
        offset += path[b].second;

        out << "\t" << path[a].first->getId() << " -> " << path[b].first->getId() <<
            " [label=" << offset << ", color=" << color << ", penwidth=" << 5 << "];" << endl;
        a = b;
        b++;
        offset = 0;
    }


}

void GraphVisualizer::visualizeGraph(Graph *G, vector<Contig*> & contigs) {

    cerr << "Creating background graph" << endl;

    FAU fau(G->size());

    for( int i=0; i<G->size(); i++ ){ // here i make union of all connected components
        int a = i;
        for( int k=0; k<(*G)[i].size(); k++ ){
            int b = (*G)[i][k].first;

            auto path = G->getContractedEdgePath( a,b );

            int x = a;
            int y = -1;
            for( auto p : path ){
                y = p.first;
                fau.Union( x, y );
                x = y;
            }
        }
    }

    cerr << "fau created, creating components to consider" << endl;

    set<int> componentsToConsider;
    for( auto ctg : contigs ){
        int compId = fau.Find( ctg->getContainedReads().back().first->getId() );
        componentsToConsider.insert( compId );
    }

    cerr << "your components to consider size(): " << componentsToConsider.size() << endl;
//    for(auto p : componentsToConsider) cerr << p << " "; cerr << endl;

    cerr << "components to consider created, writing graph background" << endl;









    for( int i=0; i<G->size(); i++ ){
        if( componentsToConsider.count( fau.Find(i) ) ) {

            for( int k=0; k<(*G)[i].size(); k++ ){
                int a = i;
                int b = (*G)[i][k].first;

                if( relevantNodes.count(a) == 0 && relevantNodes.count(b) == 0 ) continue;

                int offset = (*G)[i][k].second;
                out << "\t" << a << " -> " << b << " [label=" << offset << ", color=black, penwidth=1];" << endl;

            }
        }
    }


//    for( int i=0; i<G->size(); i++ ){ // HERE I write graph (small black lines, to show 'background' of contigs
//        int a = i;
//
////    for( int a : relevantNodes ){
//        if( componentsToConsider.count( fau.Find(a) ) ) {
//
//            if(relevantNodes.count(a)==0) continue;
//
//            for( int k=0; k<(*G)[a].size(); k++ ){
//                int b = (*G)[a][k].first;
//
//
//                auto P = G->getContractedEdgePath(a,b);
//
//                if( P.empty() ){
////                    out << "\t" << a << " -> " << b <<
////                        "[label=" << V[a][k].second << ", color=black, penwidth=1];" << endl;
//                    continue;
//                }
//
//                VPII path( P.begin(), P.end() );
//
//                int offset = 0;
//                int x = a;
//                int y = 0;
//
//                while( y < path.size() ){
//                    while( y < path.size() && relevantNodes.count( path.at(y).first ) == 0 ){
//                        offset += path[y].second;
//                        y++;
//                    }
//                    if( y == path.size() ){
//                        if( x == path.back().first ) break;
//                        else y--;
//                    }
//
//                    offset += path.at(y).second;
//
//                    out << "\t" << x << " -> " << path[y].first <<
//                        "[label=" << offset << ", color=black, penwidth=1];" << endl;
//
//                    x = path[y].first;
//                    y++;
//                    offset = 0;
//                }
//            }
//        }
//
//    }


    cerr << "Background graph created" << endl << endl;
}

void GraphVisualizer::createRelevantNodes(Graph *G, vector<Contig *> &contigs) {
    cerr << "creating relevant nodes" << endl;

    set<int> nodesInContigs;
    for(auto ctg : contigs){ // here i add the beginning and the end of the contig
        relevantNodes.insert( ctg->getContainedReads()[0].first->getId() );
        relevantNodes.insert( ctg->getContainedReads().back().first->getId() );

        for( auto p : ctg->getContainedReads() ){
            nodesInContigs.insert( p.first->getId() );
        }
    }

    Graph GRev = G->getReverseGraph();

    for( int i=0; i<G->size(); i++ ) {
        if( nodesInContigs.count(i) == 0 ) continue;

        if( (*G)[i].size() == 0 && GRev[i].size() > 0 ) relevantNodes.insert(i);
        if( (*G)[i].size() > 0 && GRev[i].size() == 0 ) relevantNodes.insert(i);
        if( (*G)[i].size() >= 2  || GRev[i].size() >= 2 ) relevantNodes.insert(i);
    }

    cerr << "There are " << relevantNodes.size() << " relevant nodes:" << endl;
//    for(auto p : relevantNodes) cerr << p << " "; cerr << endl;

    cerr << "relevant nodes created" << endl << endl;
}

void GraphVisualizer::writeWholeGraph(Graph *G, vector<Read *> &reads, string filename) {
//    out.open( Params::outStreamFileName + "_mopp" + to_string(Params::MAX_OFFSET_PARALLEL_PATHS)
//             + "_modb" + to_string(Params::MAX_OFFSET_DANGLING_BRANCHES) + "_rsoe" + to_string(Params::REMOVE_SMALL_OVERLAP_EDGES_MIN_OVERLAP) + "-" +
//             to_string(Params::REMOVE_SMALL_OVERLAP_EDGES_NUMBER_TO_RETAIN) + ".gv");

    out.open(filename + ".gv");
    assert(out.is_open());

    out << "digraph G{" << endl;

    for( int i=0; i<G->size(); i++ ){
        for( PII neigh : (*G)[i] ){
            int d = neigh.first;
            int offset = neigh.second;
            int overlap = reads[i]->size() - offset;

            string color = colors[i%C];
            out << "\t" << i << " -> " << d << " [label=\"(" << offset << "," << overlap << ")\", color=" << color << ", penwidth=" << 4 << "];" << endl;

        }
    }

    out << "}" << endl;
    out.close();

    cerr << "writing in graphviz format finished" << endl;

}



