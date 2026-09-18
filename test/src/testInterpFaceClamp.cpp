// Regression gate for the element-face coordinate clamp in
// ot::da::lagrangeInterpElementToCoord (include/daUtils.tcc).
//
// The clamp snaps a requested coordinate lying within 1e-6 of an element face
// onto that face, to absorb roundoff in the caller's coordinates. For the upper
// z face it snapped to pt_min.z() rather than pt_max.z(), so a point just below
// an element's top in z was evaluated a full element lower down. The x and y
// faces were correct; only z was wrong, and the same typo appeared in both
// lagrangeInterpElementToCoord and linear_lagrange.
//
// The probe field is linear. Lagrange interpolation of any order >= 1 reproduces
// a linear field exactly, and does so from *any* element containing the point,
// so this test does not depend on which element an ambiguous near-face point is
// assigned to -- the only way to miss is to evaluate at the wrong coordinate,
// which is precisely the defect. A wrong-element read of a smooth nonlinear
// field would be hard to distinguish from interpolation error; here it is not.
//
// Every local element is probed just inside all six of its faces. Under the bug
// the upper-z probes come back short by cz * (pt_max.z() - pt_min.z()), i.e. one
// full element height, while the other five faces pass -- so the lower-face and
// x/y-face probes double as the check that a fix does not disturb the paths that
// were already correct.
//
// Exit code 0 = pass.

#include <mpi.h>

#include <cmath>
#include <cstdio>
#include <vector>

#include "TreeNode.h"
#include "daUtils.h"
#include "dendro.h"
#include "mesh.h"
#include "octUtils.h"
#include "point.h"

// linear probe field, in domain coordinates
static const double C0 = 1.0, CX = 2.0, CY = 3.0, CZ = 5.0;

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;

    int rank, npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    m_uiMaxDepth                 = 8;
    const unsigned int regLev    = 3;
    const unsigned int eOrder    = 6;

    const double d_min           = -0.5;
    const double d_max           = 0.5;

    // how far inside a face to probe. must be under the routine's 1e-6 clamp
    // window so the clamp actually fires, and far under the element size so the
    // point is unambiguously interior.
    const double eps             = 5e-7;

    // Two scales matter here, and the tolerance has to sit between them.
    //
    // Legitimate: the clamp snaps the probe onto the face it is near, which
    // shifts the evaluation point by at most eps and the value by at most
    // max|c|*eps = 2.5e-6. That is the routine working as intended.
    //
    // The defect: snapping onto a *different* face moves the evaluation point a
    // full element, changing the value by |c|*h. At regLev 3 over a unit domain
    // h = 0.125, so the z error is 0.625 -- roughly 250000x the legitimate
    // shift. Anything in between would be a new and different bug.
    const double tol             = 1e-4;

    _InitializeHcurve(m_uiDim);

    const Point dom_min(d_min, d_min, d_min);
    const Point dom_max(d_max, d_max, d_max);

    const double gridRange   = (double)(1u << m_uiMaxDepth);
    const double domainRange = d_max - d_min;
    const double g2d         = domainRange / gridRange;

    std::vector<ot::TreeNode> balOct;
    createRegularOctree(balOct, regLev, m_uiDim, m_uiMaxDepth, comm);

    ot::Mesh* mesh = new ot::Mesh(balOct, 1, eOrder, comm);
    mesh->setDomainBounds(dom_min, dom_max);

    // sampled in octree coordinates, evaluated in domain coordinates
    std::function<double(double, double, double)> func =
        [&](const double x, const double y, const double z) {
            const double xx = (x / gridRange) * domainRange + d_min;
            const double yy = (y / gridRange) * domainRange + d_min;
            const double zz = (z / gridRange) * domainRange + d_min;
            return C0 + CX * xx + CY * yy + CZ * zz;
        };

    std::vector<double> funcVal;
    mesh->createVector(funcVal, func);
    mesh->performGhostExchange(funcVal);

    long n_bad = 0, n_checked = 0;

    if (mesh->isActive()) {
        const std::vector<ot::TreeNode>& elements = mesh->getAllElements();

        static const char* face_name[6] = {"x-lo", "x-hi", "y-lo",
                                           "y-hi", "z-lo", "z-hi"};

        std::vector<double> coords;
        std::vector<int> face_of;  // which face each probe came from

        for (unsigned int ele = mesh->getElementLocalBegin();
             ele < mesh->getElementLocalEnd(); ele++) {
            const ot::TreeNode& oct = elements[ele];

            const double exmin = d_min + oct.minX() * g2d;
            const double eymin = d_min + oct.minY() * g2d;
            const double ezmin = d_min + oct.minZ() * g2d;
            const double exmax = d_min + oct.maxX() * g2d;
            const double eymax = d_min + oct.maxY() * g2d;
            const double ezmax = d_min + oct.maxZ() * g2d;

            const double cx = 0.5 * (exmin + exmax);
            const double cy = 0.5 * (eymin + eymax);
            const double cz = 0.5 * (ezmin + ezmax);

            // six probes, each just inside one face, centered in the other two
            const double px[6] = {exmin + eps, exmax - eps, cx, cx, cx, cx};
            const double py[6] = {cy, cy, eymin + eps, eymax - eps, cy, cy};
            const double pz[6] = {cz, cz, cz, cz, ezmin + eps, ezmax - eps};

            for (int f = 0; f < 6; f++) {
                coords.push_back(px[f]);
                coords.push_back(py[f]);
                coords.push_back(pz[f]);
                face_of.push_back(f);
            }
        }

        const unsigned int nPts = face_of.size();
        std::vector<double> out(nPts, 0.0);
        std::vector<unsigned int> validIndex;

        Point grid_limit[2]   = {Point(0.0, 0.0, 0.0),
                                 Point(gridRange, gridRange, gridRange)};
        Point domain_limit[2] = {dom_min, dom_max};

        ot::da::interpolateToCoords(mesh, funcVal.data(), coords.data(),
                                    coords.size(), grid_limit, domain_limit,
                                    out.data(), validIndex);

        double worst[6] = {0, 0, 0, 0, 0, 0};
        long reported   = 0;
        for (unsigned int i = 0; i < validIndex.size(); i++) {
            const unsigned int p = validIndex[i];
            const double x = coords[3 * p + 0];
            const double y = coords[3 * p + 1];
            const double z = coords[3 * p + 2];

            const double expect = C0 + CX * x + CY * y + CZ * z;
            const double diff   = fabs(out[p] - expect);

            n_checked++;
            if (diff > worst[face_of[p]]) worst[face_of[p]] = diff;
            if (diff > tol) {
                n_bad++;
                if (reported < 5) {  // keep the log readable on failure
                    printf(
                        "  [rank %d] %s face: probe (% .9f, % .9f, % .9f) "
                        "expected % .9f got % .9f  diff %.3e\n",
                        rank, face_name[face_of[p]], x, y, z, expect, out[p],
                        diff);
                    reported++;
                }
            }
        }

        if (rank == 0) {
            printf("  worst deviation per face (tolerance %.1e):\n", tol);
            for (int f = 0; f < 6; f++)
                printf("    %-5s %.3e\n", face_name[f], worst[f]);
        }
    }

    long n_bad_g = 0, n_checked_g = 0;
    MPI_Reduce(&n_bad, &n_bad_g, 1, MPI_LONG, MPI_SUM, 0, comm);
    MPI_Reduce(&n_checked, &n_checked_g, 1, MPI_LONG, MPI_SUM, 0, comm);

    int status = 0;
    if (!rank) {
        printf("testInterpFaceClamp: checked %ld near-face probes, %ld wrong\n",
               n_checked_g, n_bad_g);
        if (n_checked_g == 0) {
            printf("FAIL: no probes were interpolated -- test is not testing "
                   "anything\n");
            status = 1;
        } else if (n_bad_g != 0) {
            printf("FAIL: element-face coordinate clamp is evaluating at the "
                   "wrong coordinate\n");
            status = 1;
        } else {
            printf("PASS\n");
        }
    }
    MPI_Bcast(&status, 1, MPI_INT, 0, comm);

    delete mesh;
    MPI_Finalize();
    return status;
}
