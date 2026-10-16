// Regrouping the blocks of an existing mesh (what ENUTS does after a remesh,
// with element time levels as tags) must leave the element-to-block unzip map
// consistent with the new block list, and unzip must stay correct.
//
// Usage: testBlockRegroupE2B [maxDepth=6] [eleOrder=6]
// Exit code 0 = pass.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <vector>

#include "TreeNode.h"
#include "dendro.h"
#include "mesh.h"
#include "meshUtils.h"
#include "mpi.h"

static const double D_MIN = -1.0;
static const double D_MAX = 1.0;

static double poly(double x, double y, double z) {
    return 1.0 + x - 2.0 * y + 0.5 * z + x * y - z * z;
}

// local failures, printed once per rank
static unsigned int check_e2b(const ot::Mesh* mesh) {
    const std::vector<ot::Block>& blkList   = mesh->getLocalBlockList();
    const std::vector<unsigned int>& map    = mesh->getE2BUnzipMap();
    const std::vector<unsigned int>& counts = mesh->getE2BUnzipMapCounts();
    const std::vector<unsigned int>& offset = mesh->getE2BUnzipMapOffsets();
    const size_t nAll                       = mesh->getAllElements().size();

    if (counts.size() != nAll || offset.size() != nAll) return 1;

    std::vector<std::vector<unsigned int>> want(nAll);
    std::vector<unsigned int> eid;
    for (unsigned int b = 0; b < blkList.size(); b++) {
        mesh->blkUnzipElementIDs(b, eid);
        for (unsigned int e = blkList[b].getLocalElementBegin();
             e < blkList[b].getLocalElementEnd(); e++)
            want[e].push_back(b);
        for (unsigned int e : eid) want[e].push_back(b);
    }

    unsigned int bad = 0;
    for (size_t e = 0; e < nAll; e++) {
        if (counts[e] != want[e].size() ||
            (size_t)offset[e] + counts[e] > map.size()) {
            bad++;
            continue;
        }
        std::vector<unsigned int> got(map.begin() + offset[e],
                                      map.begin() + offset[e] + counts[e]);
        std::sort(got.begin(), got.end());
        std::sort(want[e].begin(), want[e].end());
        if (got != want[e]) bad++;
    }
    return bad;
}

// unzip points of non-boundary blocks that disagree with poly()
static unsigned int check_unzip(ot::Mesh* mesh, const double* cg) {
    std::vector<double> uz(mesh->getDegOfFreedomUnZip(), 0.0);
    mesh->unzip(cg, uz.data(), 1);

    const double scale = (D_MAX - D_MIN) / (double)(1u << m_uiMaxDepth);
    const std::vector<ot::Block>& blkList = mesh->getLocalBlockList();
    unsigned int bad                      = 0;
    for (unsigned int b = 0; b < blkList.size(); b++) {
        if (blkList[b].getBlkNodeFlag()) continue;
        const ot::TreeNode blkNode = blkList[b].getBlockNode();
        const unsigned int lx      = blkList[b].getAllocationSzX();
        const unsigned int ly      = blkList[b].getAllocationSzY();
        const unsigned int lz      = blkList[b].getAllocationSzZ();
        const unsigned int pW      = blkList[b].get1DPadWidth();
        const double hx            = blkList[b].computeGridDx();
        const double hy            = blkList[b].computeGridDy();
        const double hz            = blkList[b].computeGridDz();
        const DendroIntL off       = blkList[b].getOffset();
        for (unsigned int k = 0; k < lz; k++)
            for (unsigned int j = 0; j < ly; j++)
                for (unsigned int i = 0; i < lx; i++) {
                    const double x =
                        (blkNode.minX() + (double(i) - pW) * hx) * scale +
                        D_MIN;
                    const double y =
                        (blkNode.minY() + (double(j) - pW) * hy) * scale +
                        D_MIN;
                    const double z =
                        (blkNode.minZ() + (double(k) - pW) * hz) * scale +
                        D_MIN;
                    const double v = uz[off + (k * ly + j) * lx + i];
                    if (!(std::fabs(v - poly(x, y, z)) < 1e-10)) bad++;
                }
    }
    return bad;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    m_uiMaxDepth              = (argc > 1) ? std::atoi(argv[1]) : 6;
    const unsigned int eOrder = (argc > 2) ? std::atoi(argv[2]) : 6;

    _InitializeHcurve(m_uiDim);

    // refined around one off-center point, so levels (and blocks) vary
    std::function<double(double, double, double)> fr = [](double x, double y,
                                                          double z) {
        const double s  = (D_MAX - D_MIN) / (double)(1u << m_uiMaxDepth);
        const double xx = x * s + D_MIN - 0.3;
        const double yy = y * s + D_MIN - 0.2;
        const double zz = z * s + D_MIN + 0.1;
        return std::exp(-40.0 * (xx * xx + yy * yy + zz * zz));
    };
    std::vector<ot::TreeNode> tmpNodes;
    function2Octree(fr, tmpNodes, m_uiMaxDepth, 1e-3, eOrder, comm);
    ot::Mesh* mesh = ot::createMesh(tmpNodes.data(), tmpNodes.size(), eOrder,
                                    comm, 1, ot::SM_TYPE::FDM);
    mesh->setDomainBounds(Point(D_MIN, D_MIN, D_MIN),
                          Point(D_MAX, D_MAX, D_MAX));

    unsigned int fails = 0;
    auto expect        = [&](bool ok, const char* what) {
        unsigned int l = ok ? 0 : 1, g = 0;
        MPI_Allreduce(&l, &g, 1, MPI_UNSIGNED, MPI_MAX, comm);
        if (g) {
            if (!rank) std::printf("  FAIL: %s\n", what);
            fails++;
        }
    };

    std::function<void(double, double, double, double*)> fcg =
        [](double x, double y, double z, double* v) { v[0] = poly(x, y, z); };
    double* cg = mesh->createCGVector<double>(fcg, 1);
    if (mesh->isActive()) {
        mesh->readFromGhostBegin(cg, 1);
        mesh->readFromGhostEnd(cg, 1);
    }

    unsigned int nblk0 = 0, nblk1 = 0, nmap1 = 0, nmap2 = 0;
    unsigned int e2b_bad = 0, uz_bad = 0, e2b_bad2 = 0;
    if (mesh->isActive()) {
        e2b_bad += check_e2b(mesh);
        uz_bad += check_unzip(mesh, cg);
        nblk0                      = mesh->getLocalBlockList().size();

        // per-element tags, as ENUTS passes its element time levels;
        // alternating in x so every multi-element block has to split
        const ot::TreeNode* pNodes = mesh->getAllElements().data();
        std::vector<unsigned int> tags(mesh->getNumLocalMeshElements());
        for (unsigned int e = 0; e < tags.size(); e++) {
            const ot::TreeNode& t = pNodes[mesh->getElementLocalBegin() + e];
            tags[e] = 1u + ((t.minX() >> (m_uiMaxDepth - t.getLevel())) & 1u);
        }

        mesh->performBlocksSetup(0, tags.data(), tags.size());
        nblk1 = mesh->getLocalBlockList().size();
        nmap1 = mesh->getE2BUnzipMap().size();
        e2b_bad += check_e2b(mesh);
        uz_bad += check_unzip(mesh, cg);

        mesh->performBlocksSetup(0, tags.data(), tags.size());
        nmap2    = mesh->getE2BUnzipMap().size();
        e2b_bad2 = check_e2b(mesh);
    }

    unsigned int nblk[2] = {nblk0, nblk1}, nblk_g[2];
    MPI_Allreduce(nblk, nblk_g, 2, MPI_UNSIGNED, MPI_SUM, comm);
    if (!rank)
        std::printf(
            "testBlockRegroupE2B npes=%d maxDepth=%u eO=%u blocks %u -> "
            "%u\n",
            npes, m_uiMaxDepth, eOrder, nblk_g[0], nblk_g[1]);

    expect(nblk_g[0] != nblk_g[1], "tags did not change the block list");
    expect(e2b_bad == 0, "e2b map disagrees with the block list");
    expect(uz_bad == 0, "unzip after regroup is wrong");
    expect(e2b_bad2 == 0 && nmap2 == nmap1,
           "repeating the same regroup changed the e2b map");

    mesh->destroyVector(cg);
    delete mesh;

    if (!rank) std::printf(fails ? "FAILED (%u)\n" : "PASS\n", fails);
    MPI_Finalize();
    return fails ? 1 : 0;
}
