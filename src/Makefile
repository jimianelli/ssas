LIB:=lib/tailcompression.cpp lib/nloglike.cpp lib/spr.cpp lib/multinomial.cpp lib/robust_multi.cpp lib/dirichlet.cpp lib/equilibrium.cpp
INC:=include/*.h include/*.hpp

ifdef ComSpec
    RM=del /F /Q
    COPY=copy
else
    RM=rm -rf
    COPY=cp
endif

all: build/debug/ssas build/release/ssas

build/debug/ssas: build/debug ssas.tpl 
	@admb -g ssas.tpl 
	@$(COPY) ssas build/debug/
	@$(RM) ssas

build/release/ssas: build/release ssas.tpl 
	@admb -g ssas.tpl 
	@$(COPY) ssas build/release/
	@$(RM) ssas

build/debug:
	@echo creating debug directory
	@mkdir -p build/debug

build/release:
	@echo creating release directory
	@mkdir -p build/release

clean:
	@$(RM) build
	@$(RM) ssas.cpp
	@$(RM) ssas.htp
	@$(RM) ssas.obj
