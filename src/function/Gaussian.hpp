#ifndef __OXSX_FITTABLE_GAUSSIAN__
#define __OXSX_FITTABLE_GAUSSIAN__

class FittableGaussian : public Gaussian, public FitComponent{

public:
    void
    MakeFittable(){
        fParameterManager.Clear();
        fParameterManager.AddContainer<std::vector<double> >(fMeans,
                                                             "Gaussian Means");
        fParameterManager.AddContainer<std::vector<double> >(fStdDevs,
                                                             "Gaussian st.devs");
    }

    std::vector<double>
    GetParameters() const{
        return fParameterManager.GetParameters();
    }

    std::vector<std::string>
    GetParameterNames() const{
        return fParameterManager.GetParameterNames();
    }

    void
    SetParameters(const std::vector<double>& params_){
        try{
            fParameterManager.SetParameters(params_);
        }
        catch(const ParameterCountError&){
            throw ParameterCountError("Fittable<Gaussian>",
                                      GetParameterCount(), params_.size(),
                                      Formatter() << GetNDims()
                                      << " means and sigmas"
                                      );
        }
    }
    
    size_t
    GetParameterCount() const{
        return fParameterManager.GetParameterCount();
    }

private:
    ParameterManager fParameterManager;
};


#endif
