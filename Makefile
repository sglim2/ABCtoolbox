subdirs = source/ABCestimator \
          source/ABCsampler \
          source/cumuldens \
          source/glm \
          source/linearTransformer \
          source/strStats

.PHONY: $(subdirs)

all: $(subdirs)
clean: $(subdirs)

$(subdirs):
	$(MAKE) -C $@ $(MAKECMDGOALS)
	

