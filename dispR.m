function res = dispR(om,acoef)
global Gamma K2 k bbeta cr2 Nr2 f;
res = om.*(cr2*(acoef.*acoef+Gamma*Gamma) + Nr2 - om.*om).*(f*f - (om+k*bbeta/K2).^2) ...
      - cr2*K2*(om+k*bbeta/K2).*(om.*om-Nr2);
end