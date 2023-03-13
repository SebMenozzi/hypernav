#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>

void __attribute__((import_module("uniforms"))) __attribute__((import_name("moebiusTransformation"))) moebiusTransformation(float fv[static 4]);
void __attribute__((import_module("uniforms"))) __attribute__((import_name("distanceBetweenLevels"))) setDistanceBetweenLevels(float f);
void __attribute__((import_module("uniforms"))) __attribute__((import_name("cellHeight"))) setCellHeight(uint16_t node);
void __attribute__((import_module("uniforms"))) __attribute__((import_name("childrenWidth"))) setChildrenWidth(float width);
void __attribute__((import_module("uniforms"))) __attribute__((import_name("bulge"))) setBulge(float bulge);
void __attribute__((import_module("uniforms"))) __attribute__((import_name("baseNode"))) setBaseNode(uint16_t node);

void __attribute__((import_module("textures"))) __attribute__((import_name("treeDataSubImage"))) treeDataSubImage(int x, int y, int w, int h, uint16_t *data);

void __attribute__((import_module("env"))) __attribute__((import_name("scheduleTick"))) scheduleTick(void);
void __attribute__((import_module("env"))) __attribute__((import_name("consoleLog"))) consoleLog(double);

struct moebius {
    complex double a;
    complex double b;
};

struct treeData {
    uint16_t parent;
    uint16_t index;
    uint16_t prev;
    uint16_t next;
    uint16_t children[4];
    uint16_t padding[16];
};

struct deceleration {
    struct moebius current;
    struct moebius target;
    double factor;
};

enum panType {
    PanTypeNone,
    PanTypeScroll,
    PanTypeSwipeBack,
    PanTypeSwipeForward,
};

static double distanceBetweenLevels;
static double cellHeight;
static double childrenWidth;
static double bulge;

static const struct moebius identity = { 1, 0 };

static struct moebius base = { 1, 0 };
static uint16_t baseNode;

static double scrollAmount[256];
static uint16_t scrollNode[256];
static int focusedScroll;
static double scrollTarget;
static double swipeBackScrollLimit;

static complex double panStart;
static complex double panPrevious;
static complex double panPrevious2;
static complex double panCurrent;
static double durationBetweenPanPreviousAndCurrent;
static double durationBetweenPanPrevious2AndCurrent;
static double durationSincePan;

static double mouseDownX;
static double mouseDownY;

static bool willPan;
static bool interactionIsTouch;
static bool interactionStoppedScrolling;
static enum panType panType;

static struct deceleration lateralDeceleration = { { 1, 0 }, { 1, 0 }, 0.985 };

static struct treeData treeData[4096];
static uint16_t numberOfNodes;

static uint16_t focusedNode;

static double viewportWidth;
static double viewportHeight;

static bool moebiusEqual(struct moebius m1, struct moebius m2)
{
    return m1.a == m2.a && m1.b == m2.b;
}
static struct moebius compose(struct moebius m1, struct moebius m2)
{
    complex double c = m2.a + m1.b * conj(m2.b);
    return (struct moebius){
        .a = m1.a * c / (m2.a * m2.b * conj(m1.b) + 1),
        .b = (m2.a * m2.b + m1.b) / c
    };
}
static struct moebius inverse(struct moebius m)
{
    return (struct moebius){ .a = conj(m.a), .b = -m.b * m.a };
}
static complex double apply(struct moebius m, complex double z)
{
    return m.a * (z + m.b) / (conj(m.b) * z + 1);
}
static struct moebius originTranslation(complex double z)
{
    return (struct moebius){ .a = 1, .b = z };
}
static struct moebius mpow(struct moebius m, double k)
{
    if (cabs(m.a - 1) < 1e-6) {
        // xx unify these cases somehow?
        double len = cabs(m.b);
        if (len < 1e-6)
            return identity;
        double Xk = pow(1 - len, k);
        double Yk = pow(1 + len, k);
        return (struct moebius){
            .a = 1,
            .b = m.b / len * (Yk - Xk) / (Xk + Yk)
        };
    }
    complex double d = csqrt(4 * m.a * m.b * conj(m.b) + (m.a - 1) * (m.a - 1));
    complex double A = d - m.a + 1;
    complex double B = d + m.a - 1;
    complex double Xk = cpow(-d + m.a + 1, k);
    complex double Yk = cpow(d + m.a + 1, k);
    return (struct moebius){
        .a = (A * Xk + B * Yk) / (B * Xk + A * Yk),
        .b = m.a * m.b * 2 * (Yk - Xk) / (A * Xk + B * Yk)
    };
}
static complex double diskToBand(complex double disk)
{
    return catanh(disk * I) / (I * 0.785) * viewportWidth;
}
static complex double bandToDisk(complex double band)
{
    return ctanh(band / viewportWidth * (I * 0.785)) / I;
}

static int rowInBand(complex double band)
{
    return (int)floor((cimag(band) + viewportHeight * 0.5) / cellHeight);
}

static double scrollInNode(double scroll, uint16_t *node)
{
    double offset = 0;
    for (int i = 0; i < 32; ++i) {
        struct treeData data = treeData[*node];
        if (scroll + offset < 0) {
            if (data.prev == 0)
                break;
            *node = data.prev;
            offset += cellHeight * 4;
        } else if (scroll + offset >= cellHeight * 4) {
            if (data.next == 0)
                break;
            *node = data.next;
            offset -= cellHeight * 4;
        } else
            break;
    }
    return offset;
}

static struct moebius moveFromParent(double row, double rowInParent)
{
    return compose(
        originTranslation(bandToDisk(I * (cellHeight * (rowInParent + 0.5) - viewportHeight * 0.5))),
        compose(
            originTranslation(distanceBetweenLevels),
            originTranslation(bandToDisk(I * cellHeight * row))
        )
    );
}

// assumes that parent nodes always have a lower node number than their children
// and that earlier nodes in a list always have a lower node number than later
// nodes.
static struct moebius moveBetweenNodes(uint16_t nodeFrom, uint16_t nodeTo)
{
    struct moebius moebiusFrom = identity;
    struct moebius moebiusTo = identity;
    while (nodeFrom != nodeTo) {
        if (treeData[nodeTo].parent == treeData[nodeFrom].parent) {
            moebiusTo = compose(originTranslation(bandToDisk(I * cellHeight * (treeData[nodeTo].index - treeData[nodeFrom].index))), moebiusTo);
            nodeTo = nodeFrom;
        } else if (nodeFrom < nodeTo) {
            moebiusTo = compose(moveFromParent(treeData[nodeTo].index, treeData[nodeTo].parent % 4), moebiusTo);
            nodeTo = treeData[nodeTo].parent / 4;
        } else {
            moebiusFrom = compose(moveFromParent(treeData[nodeFrom].index, treeData[nodeFrom].parent % 4), moebiusFrom);
            nodeFrom = treeData[nodeFrom].parent / 4;
        }
    }
    return compose(inverse(moebiusFrom), moebiusTo);
}

static struct moebius domainTransformationStep(complex double z, uint16_t *node)
{
    struct treeData data = treeData[*node];
    uint16_t parent = data.parent / 4;
    
    complex double band = diskToBand(z);
    complex double parentBand = diskToBand(apply(originTranslation(distanceBetweenLevels), bandToDisk(band + I * cellHeight * data.index))) + I * 0.5 * cellHeight;
    if (creal(parentBand) < viewportWidth - childrenWidth || fabs(cimag(parentBand) - 0.5 * cellHeight) > 0.5 * cellHeight) {
        if (parent == 0)
            return identity;
        *node = parent;
        return moveFromParent(data.index, data.parent % 4);
    } else {
        double scrollOffset = scrollInNode(cimag(band) + viewportHeight * 0.5, node);
        if (scrollOffset != 0)
            return originTranslation(bandToDisk(I * scrollOffset));
        int row = rowInBand(band);
        band -= I * (row * cellHeight - viewportHeight * 0.5);
        uint16_t next = 0;
        if (row >= 0 && row < 4)
            next = data.children[row];
        if (creal(band) < viewportWidth - childrenWidth || next == 0)
            return identity;
        *node = next;
        return inverse(moveFromParent(0, row));
    }
}
static struct moebius domainTransformation(complex double z, uint16_t *node)
{
    struct moebius m = identity;
    for (int i = 0; i < 32; ++i) {
        struct moebius d = domainTransformationStep(apply(m, z), node);
        if (moebiusEqual(d, identity))
            break;
        m = compose(d, m);
    }
    return m;
}

static double rubberband(double scrollAmount, double constrained, double factor)
{
    double delta = constrained - scrollAmount;
    if (delta == 0)
        return scrollAmount;
    return constrained - copysign((1 - (1 / ((fabs(delta) * factor / viewportHeight) + 1))) * viewportHeight, delta);
}

static double unrubberband(double rubberbandedScrollAmount, double constrained, double factor)
{
    double delta = constrained - rubberbandedScrollAmount;
    if (delta == 0)
        return rubberbandedScrollAmount;
    return constrained - copysign(viewportHeight / factor * (1 / (1 - fabs(delta) / viewportHeight) - 1), delta);
}

static double falloff(double x)
{
    // return tanh(x);
    return 1 - 1 / exp(x);
}

static complex double panLocationForXY(double x, double y, bool hasParent)
{
    complex double band = diskToBand((x + I * y) / viewportWidth);
    if (!hasParent && x > 0)
        x = rubberband(x, 0, 0.4);
    const double ramp = 0.25;
    if (x / viewportWidth >= 1 - ramp)
        x = viewportWidth * (1 + ramp * (falloff((x / viewportWidth - 1 + ramp) / ramp) - 1));
    x *= 0.9;
    return bandToDisk(x + I * y);
    // return bandToDisk(0.2 + I * y);
}

static struct moebius moebiusForPan(complex double start, complex double current)
{
    struct moebius m = originTranslation(-start);
    complex double delta = apply(inverse(m), -current);
    return compose(inverse(m), compose(originTranslation(delta), m));
}

static void setMoebiusTransformation(struct moebius m)
{
    moebiusTransformation((float[]){ creal(m.a), cimag(m.a), creal(m.b), cimag(m.b) });
}

static bool decelerationActive(struct deceleration d)
{
    return (!moebiusEqual(d.target, identity) || !moebiusEqual(d.current, identity));
}

// returns false if nothing changed (so any dependent updates can be skipped).
static bool decelerate(struct deceleration *d, double dt)
{
    if (!decelerationActive(*d))
        return false;
    double k = pow(d->factor, 1000 * dt);
    d->current = compose(mpow(compose(d->current, inverse(d->target)), k), d->target);
    d->target = mpow(d->target, k);
    d->current.a /= cabs(d->current.a);
    d->target.a /= cabs(d->target.a);
    if (cabs(d->current.b) > 1)
        d->current.b /= cabs(d->current.b);
    if (cabs(d->target.b) > 1)
        d->target.b /= cabs(d->target.b);
    if (decelerationActive(*d))
        scheduleTick();
    return true;
}

static struct moebius decelerationOffsetForVelocity(struct moebius dm, double dt, double factor)
{
    if (dt <= 0 || moebiusEqual(dm, identity))
        return identity;
    double k = pow(factor, dt * 1000);
    return mpow(dm, 1/(1 - k));
}

static double constrainedScrollAmount(double scrollAmount, uint16_t scrollNode)
{
    uint16_t bottomNode = scrollNode;
    double bottomOffset = scrollInNode(scrollAmount + viewportHeight, &bottomNode);
    if (scrollAmount + bottomOffset + viewportHeight > 4 * cellHeight)
        scrollAmount = 4 * cellHeight - bottomOffset - viewportHeight;
    uint16_t topNode = scrollNode;
    double topOffset = scrollInNode(scrollAmount, &topNode);
    if (scrollAmount + topOffset < 0)
        scrollAmount = -topOffset;
    return scrollAmount;
}

static void applyInputTransformation(double *x, double *y)
{
#if 0
    // enable this to handle events in the disk visualization.
    complex double z = *x + I * *y;
    z = diskToBand(z / viewportWidth);
    *x = creal(z);
    *y = cimag(z);
#endif
}
static bool childForMouse(uint16_t *outNode)
{
    complex double band = diskToBand(apply(compose(moveBetweenNodes(focusedNode, baseNode), compose(base, lateralDeceleration.current)), bandToDisk(mouseDownX + I * mouseDownY)));
    // debugNumber(1, rowInBand(band));
    uint16_t originallyFocusedNode = focusedNode;
#if 0
    // this code lets you tap back into the parent node / tap to switch to sibling nodes;
    // it's not useful if you're just looking at the current list.

    if (creal(band) < -bulge && treeData[focusedNode].parent / 4) {
        focusedNode = treeData[focusedNode].parent / 4;
        while (treeData[focusedNode].prev)
            focusedNode = treeData[focusedNode].prev;
        focusedScroll--;
        scrollTarget = scrollAmount[focusedScroll];
        band = diskToBand(apply(compose(moveBetweenNodes(focusedNode, baseNode), compose(base, lateralDeceleration.current)), bandToDisk(mouseDownX + I * mouseDownY)));
        struct moebius destination = compose(moveBetweenNodes(baseNode, scrollNode[focusedScroll]), originTranslation(bandToDisk(I * scrollAmount[focusedScroll])));
        lateralDeceleration.current = compose(inverse(destination), compose(base, lateralDeceleration.current));
        lateralDeceleration.target = compose(inverse(destination), compose(base, lateralDeceleration.target));
        base = destination;
    }
#endif
    int row = rowInBand(band);
    uint16_t node = focusedNode;
    while (row < 0 && treeData[node].prev) {
        node = treeData[node].prev;
        row += 4;
    }
    while (row >= 4 && treeData[node].next) {
        node = treeData[node].next;
        row -= 4;
    }
    if (row >= 0 && row < 4 && treeData[node].children[row] && treeData[node].children[row] != originallyFocusedNode) {
        *outNode = treeData[node].children[row];
        return true;
    }
    return false;
}

__attribute__((export_name("init"))) void init(void)
{
    focusedNode = 1;
    baseNode = 1;
    focusedScroll = 1;

    scrollNode[0] = 1;
    scrollNode[1] = 1;
    
    setBaseNode(baseNode);

    struct treeWorklistItem {
        uint16_t node;
        uint16_t indexInNode;
        uint16_t depth;
    } worklist[4096] = {
        { 
            .node = 0, 
            .indexInNode = 0, 
            .depth = 1 
        }
    };

    int worklistLength = 1;
    uint16_t nextIndex = 1;
    unsigned numberOfLevels = 4;

    while (worklistLength > 0) {
        if (nextIndex > 4050)
            break;
        
        struct treeWorklistItem item = worklist[--worklistLength];
        treeData[item.node].children[item.indexInNode] = nextIndex;

        uint16_t n = item.depth;

        for (uint16_t i = 0; i < n; ++i) {
            treeData[nextIndex].parent = item.node * 4 + item.indexInNode;
            treeData[nextIndex].index = i * 4;
            treeData[nextIndex].prev = i == 0 ? 0 : nextIndex - 1;
            treeData[nextIndex].next = i + 1 == n ? 0 : nextIndex + 1;

            if (item.depth < numberOfLevels) {
                for (uint16_t j = 0; j < 4; ++j) {
                    worklist[worklistLength++] = (struct treeWorklistItem){
                        .node = nextIndex,
                        .indexInNode = j,
                        .depth = item.depth + 1
                    };
                }
            }
            
            nextIndex++;
        }
    }
    numberOfNodes = nextIndex;
    treeDataSubImage(0, 0, 6, nextIndex, (uint16_t *)treeData);
}
__attribute__((export_name("resize"))) void resize(double width, double height)
{
    viewportWidth = width;
    viewportHeight = height;
    swipeBackScrollLimit = viewportWidth * 5;
    cellHeight = 70;
    setCellHeight(cellHeight);
    childrenWidth = cellHeight * 0.5 + 40;
    distanceBetweenLevels = bandToDisk(viewportWidth - cellHeight * 0.5 + 12);
    bulge = viewportWidth / 4;
    setChildrenWidth(childrenWidth);
    setBulge(bulge);
    setDistanceBetweenLevels(distanceBetweenLevels);
}
__attribute__((export_name("setUsingTouchInput"))) void setUsingTouchInput(bool flag)
{
    interactionIsTouch = flag;
}
__attribute__((export_name("mousedown"))) void mousedown(double x, double y)
{
    applyInputTransformation(&x, &y);
    willPan = true;
    interactionStoppedScrolling = fabs(scrollTarget - scrollAmount[focusedScroll]) > 0.5;
    scrollTarget = round(scrollAmount[focusedScroll]);
    scrollAmount[focusedScroll] = scrollTarget;
    mouseDownX = x;
    mouseDownY = y;
}
__attribute__((export_name("mousemove"))) void mousemove(double x, double y)
{
    applyInputTransformation(&x, &y);
    complex double delta = (x - mouseDownX) + I * (y - mouseDownY);

    if (willPan && cabs(delta) > (interactionIsTouch ? 10 : 5)) {
        double direction = fabs(carg(delta)) / M_PI;
        // wait until the first movement to pan in order to support iOS, which
        // itself waits until a finger moves a certain distance before starting
        // to send move events.
        if (direction < 0.85) {
            if (direction < 0.25) {
                base = compose(base, lateralDeceleration.current);
                lateralDeceleration.current = identity;
                lateralDeceleration.target = identity;
                panType = PanTypeSwipeBack;
                panCurrent = panLocationForXY(x, y, treeData[focusedNode].parent);
                // this negative sign makes swiping back look nicer.  it's here
                // because i messed up the code earlier and naturally got this
                // behavior, which i liked.  at the time of writing, i'm just
                // trying to finish this, so i restored the old messed-up
                // behavior rather than trying to come up with something more
                // principled.
                panStart = -panCurrent;
                // limit the magnitude of panStart to avoid precision issues.
                complex double limit = -apply(originTranslation(bandToDisk(I * swipeBackScrollLimit)), -panStart);
                panStart = -apply(compose(moveBetweenNodes(focusedNode, scrollNode[focusedScroll]), originTranslation(bandToDisk(I * scrollAmount[focusedScroll]))), -panStart);
                if (!(cabs(limit) > cabs(panStart)))
                    panStart = limit;
            } else {
                panType = PanTypeScroll;
                double constrained = constrainedScrollAmount(scrollAmount[focusedScroll], scrollNode[focusedScroll]);
                panStart = y + unrubberband(scrollAmount[focusedScroll], constrained, 0.55);
                panCurrent = y;
            }
        } else {
            base = compose(base, lateralDeceleration.current);
            lateralDeceleration.current = identity;
            lateralDeceleration.target = identity;
            panType = PanTypeSwipeForward;
            panStart = panLocationForXY(x, y, true);
            panCurrent = panStart;
        }
        willPan = false;
        panPrevious2 = panCurrent;
        panPrevious = panCurrent;
        durationSincePan = 0;
        scheduleTick();
    }
    if (panType != PanTypeNone) {
        if (durationSincePan > 0) {
            durationBetweenPanPrevious2AndCurrent = durationBetweenPanPreviousAndCurrent + durationSincePan;
            durationBetweenPanPreviousAndCurrent = durationSincePan;
            panPrevious2 = panPrevious;
            panPrevious = panCurrent;
            durationSincePan = 0;
            scheduleTick();
        }
        if (panType == PanTypeScroll) {
            panCurrent = y;
            scrollAmount[focusedScroll] = panStart - panCurrent;

            // constrain scroll offset and apply rubberbanding.
            double constrained = constrainedScrollAmount(scrollAmount[focusedScroll], scrollNode[focusedScroll]);
            double constraintDelta = constrained - scrollAmount[focusedScroll];
            if (scrollAmount[focusedScroll] != constrained) {
                double r = rubberband(scrollAmount[focusedScroll], constrained, 0.55) - scrollAmount[focusedScroll];
                panCurrent -= r;
                scrollAmount[focusedScroll] += r;
            }

            // normalize scroll offset.
            double offset = scrollInNode(scrollAmount[focusedScroll], &scrollNode[focusedScroll]);
            scrollAmount[focusedScroll] += offset;
            panStart += offset;
            scrollTarget = scrollAmount[focusedScroll];

            // set base according to scroll amount.
            base = compose(moveBetweenNodes(baseNode, scrollNode[focusedScroll]), originTranslation(bandToDisk(I * scrollAmount[focusedScroll])));
        } else {
            base = compose(base, inverse(moebiusForPan(panStart, panCurrent)));
            panCurrent = panLocationForXY(x, y, panType == PanTypeSwipeForward || treeData[focusedNode].parent);
            base = compose(base, moebiusForPan(panStart, panCurrent));
        }
        base = compose(domainTransformation(apply(base, 1e-3), &baseNode), base);
        base.a /= cabs(base.a);
        setBaseNode(baseNode);
        setMoebiusTransformation(base);
    }
}
__attribute__((export_name("mouseup"))) void mouseup(bool commit)
{
    uint16_t childNode;
    int childRow;
    if (willPan && !interactionStoppedScrolling && commit) {
        if (childForMouse(&focusedNode)) {
            focusedScroll++;
            scrollNode[focusedScroll] = focusedNode;
            scrollAmount[focusedScroll] = 0;
            scrollTarget = 0;
            struct moebius movement = compose(moveBetweenNodes(focusedNode, baseNode), base);
            lateralDeceleration.current = compose(movement, lateralDeceleration.current);
            lateralDeceleration.target = compose(movement, lateralDeceleration.target);
            base = compose(base, inverse(movement));
        }
    }
    willPan = false;
    if (panType != PanTypeNone && commit) {
        if (panType == PanTypeSwipeBack || panType == PanTypeSwipeForward) {
            struct moebius velocity = identity;
            if (durationSincePan < 0.05)
                velocity = compose(inverse(moebiusForPan(panStart, panPrevious)), moebiusForPan(panStart, panCurrent));
            if (durationBetweenPanPreviousAndCurrent > 0 && commit) {
                if (panType == PanTypeSwipeBack && treeData[focusedNode].parent / 4 && creal(apply(velocity, 0)) < 0) {
                    focusedNode = treeData[focusedNode].parent / 4;
                    while (treeData[focusedNode].prev)
                        focusedNode = treeData[focusedNode].prev;
                    focusedScroll--;
                } else if (panType == PanTypeSwipeForward && creal(apply(velocity, 0)) > 0 && childForMouse(&focusedNode)) {
                    focusedScroll++;
                    scrollNode[focusedScroll] = focusedNode;
                    scrollAmount[focusedScroll] = 0;
                }
            }
            scrollTarget = scrollAmount[focusedScroll];
            struct moebius destination = compose(moveBetweenNodes(baseNode, scrollNode[focusedScroll]), originTranslation(bandToDisk(I * scrollAmount[focusedScroll])));
            struct moebius offset = decelerationOffsetForVelocity(velocity, durationBetweenPanPreviousAndCurrent, lateralDeceleration.factor);
            lateralDeceleration.current = compose(inverse(destination), compose(base, lateralDeceleration.current));
            lateralDeceleration.target = compose(lateralDeceleration.current, offset);
            base = destination;
        } else if (panType == PanTypeScroll) {
            double dz = panPrevious2 - panCurrent;
            if (durationBetweenPanPrevious2AndCurrent > 0 && commit && durationSincePan < 0.05) {
                double k = pow(0.998, 1000 * durationBetweenPanPrevious2AndCurrent);
                scrollTarget = round(scrollAmount[focusedScroll] + dz / (1 - k));
            } else
                scrollTarget = scrollAmount[focusedScroll];
            if (fabs(scrollTarget - scrollAmount[focusedScroll]) < 100)
                scrollTarget = scrollAmount[focusedScroll];
        }
        panType = PanTypeNone;
    }
    if (decelerationActive(lateralDeceleration) || scrollTarget != scrollAmount[focusedScroll] || scrollAmount[focusedScroll] != constrainedScrollAmount(scrollAmount[focusedScroll], scrollNode[focusedScroll]))
        scheduleTick();
}
__attribute__((export_name("keydown"))) void keydown(int keycode)
{
}
__attribute__((export_name("keyup"))) void keyup(int keycode)
{
}
__attribute__((export_name("tick"))) void tick(double timestamp)
{
    static double lastTimestamp;
    double dt = timestamp - lastTimestamp;
    lastTimestamp = timestamp;
    if (dt > 1. / 30)
        dt = 1. / 30;
    if (panType != PanTypeNone) {
        durationSincePan += dt;
        if (durationSincePan < 0.05)
            scheduleTick();
    }
    struct moebius deceleratingBase = base;
    bool decelerating = false;
    if (panType == PanTypeNone && !willPan && (scrollAmount[focusedScroll] != scrollTarget || scrollAmount[focusedScroll] != constrainedScrollAmount(scrollAmount[focusedScroll], scrollNode[focusedScroll]))) {
        double k = pow(0.998, 1000 * dt);
        double lastScrollAmount = scrollAmount[focusedScroll];
        double delta = (scrollTarget - scrollAmount[focusedScroll]) * (1 - k);
        scrollAmount[focusedScroll] += delta;
        double offset = scrollInNode(scrollAmount[focusedScroll], &scrollNode[focusedScroll]);
        scrollAmount[focusedScroll] += offset;
        scrollTarget += offset;

        double constrained = constrainedScrollAmount(scrollAmount[focusedScroll], scrollNode[focusedScroll]);
        double constraintDelta = constrained - scrollAmount[focusedScroll];
        if (scrollAmount[focusedScroll] != constrained) {
            double bounceK = pow(0.99, 1000 * dt);
            scrollAmount[focusedScroll] = scrollAmount[focusedScroll] * bounceK + constrained * (1 - bounceK);
            scrollTarget = scrollTarget * bounceK + constrained * (1 - bounceK);
            if (fabs(scrollAmount[focusedScroll] - constrained) < 0.5 && fabs(scrollTarget - constrained) < 0.5) {
                scrollAmount[focusedScroll] = constrained;
                scrollTarget = constrained;
            }
        } else if (fabs(scrollAmount[focusedScroll] - scrollTarget) < 3)
            scrollTarget = scrollAmount[focusedScroll]; // xx it didn't actually reach the target... also it's not rounded any more

        base = compose(moveBetweenNodes(baseNode, scrollNode[focusedScroll]), originTranslation(bandToDisk(I * scrollAmount[focusedScroll])));
        deceleratingBase = base;
        decelerating = true;
        if (scrollAmount[focusedScroll] != scrollTarget || scrollAmount[focusedScroll] != constrained)
            scheduleTick();
    }
    if (decelerate(&lateralDeceleration, dt)) {
        deceleratingBase = compose(deceleratingBase, lateralDeceleration.current);
        decelerating = true;
    }
    if (decelerating) {
        struct moebius d = domainTransformation(apply(deceleratingBase, 1e-3), &baseNode);
        base = compose(d, base);
        deceleratingBase = compose(d, deceleratingBase);
        base.a /= cabs(base.a);
        deceleratingBase.a /= cabs(deceleratingBase.a);
        setBaseNode(baseNode);
        setMoebiusTransformation(deceleratingBase);
    }
}