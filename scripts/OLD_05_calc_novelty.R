library(tidyverse)

calculate_novelty <-
  function(selection_set_results){
    # browser()
    genome_summary <-
      selection_set_results %>%
      mutate(genome_vectors=map(selection_set, 1),
             selection_orders=map(genome_vectors, selection_orders)) %>%
      select(selection_orders) %>%
      unnest(selection_orders) %>%
      filter(selected_order !=1) %>%
      group_by(genome_name) %>%
      summarise(mean_rank=mean(selected_order),
                median_rank=median(selected_order),
                number_selections=n(),
                best_rank=min(selected_order),
                worst_rank=max(selected_order)) %>%
      arrange(median_rank) %>%
      mutate(novelty_score1=number_selections*((1/(mean_rank + median_rank))),
             novelty_score2=(number_selections + (number_selections/(mean_rank + median_rank))),
             novelty_score3=((number_selections^2 / (mean_rank + median_rank))),
             novelty_score4=((number_selections / (median_rank))),
             novelty_score5=((number_selections^2 / (median_rank)^2)),
             novelty_score6=number_selections/max(number_selections) / (median_rank/n())) %>%
      
      filter(!(number_selections == 1 & best_rank == 1)) %>%
      transmute(asm_acc=genome_name,
                median_rank,
                number_selections,
                best_rank,
                worst_rank,
                novelty_score=novelty_score4,
                log_novelty=log(novelty_score)) %>%
      arrange(desc(novelty_score)) %>%
      ungroup() %>%
      mutate(RANK=1:n())
    
    return(genome_summary)
    
  }

selection_orders <-
  function(VECTOR){
    tibble(genome_name=VECTOR,
           selected_order=1:length(VECTOR))
  }

### need to calculate derep sets first.  Dont have a script for that yet.
derep_sets <- read_rds('output/derep_sets.rds')

novelty_info <-
  calculate_novelty(derep_sets) %>%
  write_tsv('./output/novelty_ranks.tsv')


hist(novelty_info$log_novelty)
novelty_info %>% ggplot(aes(x=median_rank, y=number_selections)) + geom_point()
novelty_info %>% ggplot(aes(x=median_rank, y=RANK)) + geom_point()

novelty_info %>% ggplot(aes(x=number_selections, y=RANK)) + geom_point()




