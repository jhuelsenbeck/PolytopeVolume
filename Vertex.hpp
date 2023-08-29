#ifndef Vertex_H
#define Vertex_H

#include "Vector.hpp"



class Vertex : public Vector {

    public:
                        Vertex(void);
                        Vertex(const Vector& v);
                        Vertex(Vertex& v);
                        Vertex(mpq_class& xq, mpq_class& yq, mpq_class& zq);
        bool            operator==(const Vertex& rhs) const;
        void            clean(void);
        Vertex*         getFrom(void) { return from; }
        Vertex*         getTo(void) { return to; }
        void            setFrom(Vertex* v) { from = v; }
        void            setTo(Vertex* v) { to = v; }
    
    private:
        Vertex*         from;
        Vertex*         to;
};


#endif
