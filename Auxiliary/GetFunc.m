function [func] = GetFunc(bankOfEstimations, iRep, iLevel)
    tmp = bankOfEstimations{iRep};
    func = tmp{iLevel};
end