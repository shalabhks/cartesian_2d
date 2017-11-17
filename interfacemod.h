      interface

        subroutine calc_spectral_radii(sr)
          use data_define
          implicit none
          real          :: sr(:,:,:)
        end subroutine calc_spectral_radii

        subroutine set_cut_cell(mcell)
          use data_define
          implicit none
          logical, allocatable  :: mcell(:)
        end subroutine set_cut_cell

        subroutine merge_cells(mcell)
          implicit none
          logical, intent(in)   ::  mcell(:)
        end subroutine merge_cells
     
      end interface
