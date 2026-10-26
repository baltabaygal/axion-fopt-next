// Adapted from checks/pockets/evolve.cpp; same binary layout.
// Cell order (x,y,z) with z fastest; per cell phi, a, pi_phi, pi_a.
template<class View,class Layout>
class PocketStateView {
    std::array<View,4> fields_;
    Layout layout_;
    int n_;
public:
    PocketStateView(std::array<View,4> fields,Layout layout,int n):fields_(fields),layout_(layout),n_(n){}
    std::size_t size() const {return 4*std::size_t(n_)*n_*n_;}
    double operator[](std::size_t scalar) const {
        const auto cell=scalar/4;
        const int i=cell/(std::size_t(n_)*n_),j=(cell/n_)%n_,k=cell%n_;
        device::IdxArray<3> p;layout_.putMemoryIndexFromSpatialLocationInto(p,i,j,k);
        for(int d=0;d<3;++d)if(p[d]<layout_.getNGhosts()||p[d]>=layout_.getNGhosts()+layout_.getLocalSizes()[d])
            throw std::runtime_error("State index outside interior");
        return fields_[scalar%4](p[0],p[1],p[2]);
    }
};
static auto state(BubbleModel& m,int n) {
    auto phi=m.fldS(0_c).getFullNDHostView(),ax=m.fldS(1_c).getFullNDHostView();
    auto pphi=m.piS(0_c).getFullNDHostView(),pax=m.piS(1_c).getFullNDHostView();
    const auto layout=m.fldS(0_c).getLayout();
    PocketStateView v(std::array{phi,ax,pphi,pax},layout,n);
    for(std::size_t scalar=0;scalar<v.size();++scalar)
        if(!std::isfinite(v[scalar]))throw std::runtime_error("Nonfinite evolved field");
    return v;
}
template<class State>
static void snapshot(const State& v,const std::filesystem::path& path,int n,double box,double u) {
    std::ofstream out(path,std::ios::binary);const std::uint32_t dimension=n;
    out.write("PKST0001",8);out.write(reinterpret_cast<const char*>(&dimension),4);
    for(double d:{box,u})out.write(reinterpret_cast<const char*>(&d),8);
    std::array<double,8192> buffer;
    for(std::size_t start=0;start<v.size();start+=buffer.size()){
        const auto count=std::min(buffer.size(),v.size()-start);
        for(std::size_t j=0;j<count;++j)buffer[j]=v[start+j];
        out.write(reinterpret_cast<const char*>(buffer.data()),count*8);
    }
    if(!out)throw std::runtime_error("Snapshot write failed");
}
// Plane through the box centre, z index n/2, all four components.
template<class State>
static void slice(const State& v,const std::filesystem::path& path,int n,double box,double u) {
    std::ofstream out(path,std::ios::binary);const std::uint32_t dimension=n;
    out.write("PKSL0001",8);out.write(reinterpret_cast<const char*>(&dimension),4);
    for(double d:{box,u})out.write(reinterpret_cast<const char*>(&d),8);
    for(int i=0;i<n;++i)for(int j=0;j<n;++j){
        const auto p=4*((std::size_t(i)*n+j)*n+n/2);
        const std::array<double,4> cell={v[p],v[p+1],v[p+2],v[p+3]};
        out.write(reinterpret_cast<const char*>(cell.data()),32);
    }
    if(!out)throw std::runtime_error("Slice write failed");
}
